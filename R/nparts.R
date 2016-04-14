nparts.atmost = function(n, k, include.zero=TRUE)
{## this computes partitions::R(k,n,include.zero)
 ## i.e., the number of restricted partitions of n into at most k parts
 ## reference: http://dlmf.nist.gov/26.9
	if(!include.zero) return(Recall(n-k, k))
	
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
	as.character(recur.fctn(n,k))
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

nparts.atmost.ubound = function(n, k, upper, include.zero=TRUE)
{## this computes the number of restricted partitions of n into at most k parts
 ## subject to the constraint that each part is at most "upper"
 ## reference: https://en.wikipedia.org/wiki/Partition_(number_theory)
 #	Andrews. (1998) The Theory of Partitions. CUP. Eqns (3.2.4), (3.2.6), (3.5.3), (3.5.4)
	if(!include.zero) return(Recall(n-k, k, upper-1L))

	if (k < upper) {tmp=k; k=upper; upper=tmp} #(3.5.3)
	n = min(upper*k-n , n) #(3.5.4)
	dims=pmax(rep(1L,3), c(n,k,upper)+1L)
	pkn.table = array(NA_character_, dim=dims)
	pkn.table[ , , 1L] = '0'
	pkn.table[1L, , ] = '1'
	if(k>=1L){
		tmp = matrix('1', dims[1], dims[3])
		tmp[lower.tri(tmp)]='0'
		pkn.table[, 2L, ] = tmp
	}
	
	recur.fctn = function(n, k, upper)
	{
		## these two switches help to reduced the number of recursive calls
		if (k < upper) {tmp=k; k=upper; upper=tmp} #(3.5.3)
		n = min(upper*k-n , n) #(3.5.4)

		if(n < 0L ) return(0L) 

		if(!is.na(tmpans <- pkn.table[n+1L, k+1L, upper+1L]))
			return(as.bigz(tmpans))

		ans = if(upper>n){
			Recall(n, k, n)
		}else if (k>n) {
			Recall(n,n,upper) 
		}else {
			Recall(n, k-1L, upper) + Recall(n-k, k, upper-1L)
		}
		pkn.table[n+1L, k+1L, upper+1L] <<- as.character(ans)
		ans
	}
	as.character(recur.fctn(n,k,upper))
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


