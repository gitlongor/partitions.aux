nparts.atmost = function(n, k, include.zero=TRUE)
{## this computes partitions::R(k,n,include.zero)
 ## i.e., the number of restricted partitions of n into at most k parts 
 ## reference: http://dlmf.nist.gov/26.9
	pkn.table = matrix(NA_character_, n+1L, k+1L)
	pkn.table[, 2L] = '1'
	pkn.table[, 1L] = '0'
	pkn.table[1L, ] = '1'  # overwrite (1,1) position correctly

	recur.fctn = function(n, k)
	{
		if(!is.na(tmpans <- pkn.table[n+1L, k+1L])) return(as.bigz(tmpans))
		ans = if(n<k) Recall(n, n) else Recall(n-k, k) + Recall(n, k-1L)
		pkn.table[n+1L, k+1L] <<- as.character(ans)
		ans
	}
	as.character(recur.fctn(n,k) - if(include.zero) 0L else recur.fctn(n,k-1L))
}

nparts.exact = function(n, k, include.zero=FALSE)
{## this computes partitions::R(k,n,include.zero)
 ## i.e., the number of restricted partitions of n in to exactly k parts
 ## reference: Arndt. Matters Computational. Eqn 16.4-7
 ##			http://mathworld.wolfram.com/PartitionFunctionP.html
 ##			http://nvlpubs.nist.gov/nistpubs/jres/74B/jresv74Bn1p1_A1b.pdf
	if(include.zero) return(Recall(n+k, k))
	
	if(n<=0L || k<=0L) return('0')
	tmp = switch(k, 1L, n%/%2L, (n*n+3L)%/%12L)
	if(!is.null(tmp)) return(as.character(tmp))
	
	pkn.table = matrix(NA_character_, n+1L, k+1L)
	pkn.table[, 4L] = as.character(divq.bigz(as.bigz(seq_len(n+1L)-1L)^2+3L, 12L))
	pkn.table[, 3L] = as.character(divq.bigz(seq_len(n+1L)-1L, 2L))
	pkn.table[, 2L] = '1'
	pkn.table[, 1L] = '0'
	pkn.table[1L, ] = '0'
	diag(pkn.table) = '1'

	recur.fctn = function(n, k)
	{
		if (k>n) return(0L)
		if(!is.na(tmpans <- pkn.table[n+1L, k+1L])) return(as.bigz(tmpans))
		ans = Recall(n-1L, k-1L) + Recall(n-k, k)
		pkn.table[n+1L, k+1L] <<- as.character(ans)
		ans
	}
	as.character(recur.fctn(n,k))
}

nparts.atmost.ubound = function(n, k, upper)
{## this computes the number of restricted partitions of n into at most k parts 
 ## subject to the constraint that each part is at most "upper"
 ## reference: https://en.wikipedia.org/wiki/Partition_(number_theory)

	pkn.table = matrix(NA_character_, n+1L, k+1L)
	pkn.table[, 2L] = '1'
	pkn.table[, 1L] = '0'
	pkn.table[1L, ] = '1'  # overwrite (1,1) position correctly
	
	recur.fctn = function(n, k, upper)
	{
		if(n < 0L) n = 0L
		if(upper < 0L || (n>0L && upper * k < n)) return(0L)
		if(upper == 0L) return( n==0L )
		if(upper >= n) return(as.bigz(nparts.atmost(n, k, include=TRUE)))
		#if(!is.na(tmpans <- pkn.table[n+1L, k+1L])) return(as.bigz(tmpans))
		ans = Recall(n, k-1L, upper) + Recall(n-k, k, upper-1L)
		#pkn.table[n+1L, k+1L] <<- as.character(ans)
		ans
	}
	as.character(recur.fctn(n,k,upper))
}

nparts.atmost.ubound = function(n, k, upper)
{## this computes the number of restricted partitions of n into at most k parts 
 ## subject to the constraint that each part is at most "upper"
 ## reference: https://en.wikipedia.org/wiki/Partition_(number_theory)

	pkn.table = array(NA_character_, dim=c(n+1L, k+1L, upper+1L))
	#pkn.table[, 2L] = '1'
	#pkn.table[, 1L] = '0'
	#pkn.table[1L, ] = '1'  # overwrite (1,1) position correctly
	
	recur.fctn = function(n, k, upper)
	{
		if(n < 0L) return(0L) #?
		if(n==0L) {
			if(k==0L && upper>=0L) return(1L) else
			if(k>0) return(upper>=0L) else
			return(0L)
		}
		if(k==1L) return(upper>=n)
		if(upper < 0L || (n>0L && upper * k < n)) return(0L)
		if(upper == 0L) return( n==0L )
		if(upper >= n) return(as.bigz(nparts.atmost(n, k, include=TRUE)))
		if(!is.na(tmpans <- pkn.table[n+1L, k+1L, upper+1L])) return(as.bigz(tmpans))
		ans = if(k>n) Recall(n,n,upper) else Recall(n, k-1L, upper) + Recall(n-k, k, upper-1L)
		pkn.table[n+1L, k+1L, upper+1L] <<- as.character(ans)
		ans
	}
	as.character(recur.fctn(n,k,upper))
}

restrictedpartition=function(n, k=n, values=seq_leng(n))
{

	values=unique(as.integer(values)); 
	values=values[values>0L & values<=n]; 
	bounds=range(values)
	if(bounds[1L] != 1L){
		values=values - bounds[1L] + 1L 
		n = n - k * bounds[1L]
		max.v = diff(bounds)
	}else max.v = bounds[2L]

	n=as.integer(n); k=as.integer(k); 
	
	## nout is an upper bound on the number of partions
	nout=as.integer(nparts.atmost.ubound(n,k,max.v))

	## fgp is based on the partition-en in fxt library
	## it finds all partitions with each part in "values" (excluding 0)
	## partitions with more than k parts are ignored
	ans=.Call(fpg, n, k, values, nout)
	if(length(values)!=n) ans=ans[,.colSums(ans, k, nout)>0,drop=FALSE]
	
	if(bounds[1L] != 1L) ans=ans+bounds[1L]
	ans
}


if(FALSE){
	## pre-compute a large nparts.atmost table (0<=n<=1000; 0<=k<=1000)
	library(gmp)
	recur.fctn.char=function(n, k)
	{## a variant updating a global (char) table
		if(!is.na(tmpans <- nparts.atmost.table[n+1L, k+1L])) return(as.bigz(tmpans))
		
		if(n<k){
			ans=Recall(n, n)
		}else 
			ans = Recall(n-k, k) + Recall(n, k-1L)
		nparts.atmost.table[n+1L, k+1L] <<- as.character(ans)
		return(ans)
	}
	
	maxn=maxk=1000L
	nparts.atmost.table=matrix(NA_character_, maxn+1L, maxk+1L)
	nparts.atmost.table[, 2L]='1'
	nparts.atmost.table[, 1L]='0'
	nparts.atmost.table[1L, ]='1'  # overwrite (1,1) position correctly
	
	all.ind=as.matrix(expand.grid(seq_len(maxn+1L), seq_len(maxk+1L)))
	ord=order(rowSums(all.ind), all.ind[,1L])
	idx=all.ind[ord,]
	for(i in seq_len(nrow(idx))){
		if(!is.na(nparts.atmost.table[idx[i,1L], idx[i,2L]])) next
		recur.fctn.char(idx[i,1L]-1L, idx[i,2L]-1L) ->tmp
		cat(i,'\t', idx[i,]-1L, '\t', as.character(tmp), '\n')
	}
	save(nparts.atmost.table, file='nparts.atmost.table.char.RData',compress='xz')
}


if(FALSE){
	## pre-compute a large nparts.exact table (0<=n<=1000; 0<=k<=1000)
	library(gmp)
	recur.fctn.char=function(n, k)
	{	
		if (k>n) return(0L)
		if(!is.na(tmpans <- nparts.exact.table[n+1L, k+1L])) return(as.bigz(tmpans))
		ans = Recall(n-1L, k-1L) + Recall(n-k, k)
		nparts.exact.table[n+1L, k+1L] <<- as.character(ans)
		ans
	}
	
	maxn=maxk=1000L
	nparts.exact.table=matrix(NA_character_, maxn+1L, maxk+1L)
	nparts.exact.table[, 2L] = '1'
	nparts.exact.table[, 1L] = '0'
	nparts.exact.table[1L, ] = '0'
	diag(nparts.exact.table) = '1'
	
	idx=as.matrix(expand.grid(seq_len(maxn+1L), seq_len(maxk+1L)))
	for(i in seq_len(nrow(idx))){
		if(idx[i, 1L] < idx[i, 2L]) {
			nparts.exact.table[idx[i,1L], idx[i,2L]]='0'
			next
		}
		if(!is.na(nparts.exact.table[idx[i,1L], idx[i,2L]])) next
		recur.fctn.char(idx[i,1L]-1L, idx[i,2L]-1L) ->tmp
		cat(i,'\t', idx[i,]-1L, '\t', as.character(tmp), '\n')
	}
	save(nparts.exact.table, file='nparts.exact.table.char.RData',compress='xz')
}
