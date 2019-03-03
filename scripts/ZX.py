# python modules
##python ZX.py multi_reads.sam >multi_reads.cntTable

import math 
import sys
import array
import operator
import pysam

def main():
	samfile = pysam.AlignmentFile(sys.argv[1],'r')
	references = samfile.references
	prev_read_name =''
	multi_align=[]
	multi_reads=[]
	multi_counts_TE={}
	rlen=[]
	num=0
	for aligned_read in samfile.fetch():	
		cur_read_name = aligned_read.query_name
		if prev_read_name == "":
			if aligned_read.get_tag("NH")==1:
				multi_align.append(references[aligned_read.tid])
				multi_reads.append(multi_align)
				multi_align=[]
				rlen.append(aligned_read.query_length)
				if references[aligned_read.tid] not in multi_counts_TE:
					multi_counts_TE[references[aligned_read.tid]] =1
				else:
					multi_counts_TE[references[aligned_read.tid]] +=1
				continue
			else:
				multi_align.append(references[aligned_read.tid])
				length=aligned_read.query_length
				multi_counts_TE[references[aligned_read.tid]]=1/float(aligned_read.get_tag("NH"))
				prev_read_name = cur_read_name			
				continue
			
		elif cur_read_name == prev_read_name :
			multi_align.append(references[aligned_read.tid])
			length=aligned_read.query_length
			if references[aligned_read.tid] not in multi_counts_TE:
				multi_counts_TE[references[aligned_read.tid]] =1/float(aligned_read.get_tag("NH"))
			else:
				multi_counts_TE[references[aligned_read.tid]] +=1/float(aligned_read.get_tag("NH"))
			prev_read_name = cur_read_name			
			continue
		
		else:
			prev_read_name = cur_read_name
			multi_reads.append(multi_align)
			rlen.append(length)
			multi_align=[]
			multi_align.append(references[aligned_read.tid])
			if references[aligned_read.tid] not in multi_counts_TE:
				multi_counts_TE[references[aligned_read.tid]] =1/float(aligned_read.get_tag("NH"))
			else:
				multi_counts_TE[references[aligned_read.tid]] +=1/float(aligned_read.get_tag("NH"))
	
	multi_counts=[0.0]*len(multi_counts_TE)
	TE_trackback={}
	num_trackback={}
	for te in multi_counts_TE:
		multi_counts[num]=multi_counts_TE[te]
		TE_trackback[te]=num
		num_trackback[num]=te
		num +=1	
		
	mappability={}
	map_file = open(sys.argv[2],'r')
	for line in map_file:
		line=line.strip()
		items=line.split('\t')
		mappability[items[0]]=float(items[1])

	uniq_counts={}
	uniq_file = open(sys.argv[3],'r')
	for line in uniq_file:
		line=line.strip()
		items=line.split('\t')
		uniq_counts[items[0]]=float(items[1])	
	
	te_multi_counts =[0]*len(multi_counts_TE)
	te_multi_counts=EM(multi_reads,multi_counts,TE_trackback,rlen,mappability,uniq_counts)
	for i in range(len(te_multi_counts)):
		print num_trackback[i],"\t",int(te_multi_counts[i])

def normalizeMeans(meansIn):
	total = sum(meansIn)
	meansOut = [0]*len(meansIn)
	sys.stderr.write("total means = " + str(total) +"\n")	 
	if total > 0 :
		meansOut = map(lambda x: 1.0*x/total, meansIn)
	return meansOut
	
def M_step(meansIn,multi_reads,TE_trackback,rlen,mappability,uniq_counts):
	meansOut = [0] *len(meansIn)
	multi_counts=E_step(meansIn,multi_reads,TE_trackback,rlen,mappability,uniq_counts)
	for tid in range(len(meansIn)) :
		meansOut[tid] = multi_counts[tid]						  
	meansOut = normalizeMeans(meansOut)
	#print meansOut
	#sys.stderr.write("after normalization toatl means = "+str(sum(meansOut))+"\n") ##################
   
	return meansOut

	
def EM(multi_reads,multi_counts,TE_trackback,rlen,mappability,uniq_counts,numItr=50,OPT_TOL=0.0001):
	#initializaiton
	means0 = []
	for tid in range(len(multi_counts)) :
		means0.append(multi_counts[tid])
	means0 = normalizeMeans(means0) 
	#print means0
	#################
	cur_iter = 0
	t_size = len(multi_counts)
	r = [0]*t_size
	outerIteration = 1
	while cur_iter < numItr :
		cur_iter += 1
		
		means1 = M_step(means0,multi_reads,TE_trackback,rlen,mappability,uniq_counts)					  
		
		for tid in range(len(means0)) :
			r[tid] = means1[tid] - means0[tid]
			
		rNorm = math.sqrt(sum(map(operator.mul,r,r)))
		
		if rNorm < OPT_TOL :
			means1 = means0
			#sys.stderr.write("rNome = OPT_TOL \n")
			break
		
		means1, means0 = means0, means1
			
	if cur_iter >= numItr :
		sys.stderr.write("not converge.....\n")
	else :
		sys.stderr.write("converge at iteration " + str(cur_iter)+"\n")	   
	new_multi_counts = E_step(means0,multi_reads,TE_trackback,rlen,mappability,uniq_counts)
	#print new_multi_counts
	return new_multi_counts
			
def E_step(meansIn,multi_reads,TE_trackback,rlen,mappability,uniq_counts):
	
	multi_counts = [0] *len(meansIn)
	
	sys.stderr.write("num of multi reads = "+str(len(multi_reads))+"\n")
	
	for f in range(len(multi_reads)) :
		totalMass = 0.0		
		for te in multi_reads[f] :			
			te_id = TE_trackback[te]
			totalMass += meansIn[te_id]*uniq_counts[te]/mappability[te]  #E_step
		
		if totalMass > 0.0 :
				  norm = 1.0 / totalMass
		else :
				  norm = 0.0
				  
		for te in multi_reads[f] :
			te_id = TE_trackback[te]	
			multi_counts[te_id] += meansIn[te_id] * norm *uniq_counts[te]/mappability[te]
	
	#sys.stderr.write("total multi counts = "+ str(sum(multi_counts))+"\n")
	return multi_counts
		
if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("User interrupt !\n")
		sys.exit(0)				