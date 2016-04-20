

restrictedpartitionFXT=function(n, k=n, values=seq_len(n))
{

	values=sort.int(unique(as.integer(values)), decreasing=FALSE);
	values=values[values>0L & values<=n];
	bounds=range(values)
	max.v = bounds[2L]

	n=as.integer(n); k=as.integer(k);

	## nout is an upper bound on the number of partions
	nout=as.integer(nparts.atmost.ubound(n,k,max.v))

	## fgp is based on the partition-en in fxt library
	## it finds all partitions with each part in "values" (excluding 0)
	## partitions with more than k parts are ignored
	.Call(fpg, n, k, values, nout)
}

restrictedpartitionR=function(n, k=n, values=seq_len(n))
{ ## recursive implementation
	values=c(0L, sort(unique(values[values>0 & values<=n]), decreasing=TRUE))
	nvals=length(values)
	maxv=max(values)
	out=matrix(NA_integer_, k, as.numeric(nparts.atmost.ubound(n,k,maxv)))
	nsol=0L

	x=rep(n, nvals)
	ct=cumct=rep(0L, nvals)

	recur.fctn=function(lev)
	{
		if(lev > nvals) return(NULL)
		maxct = min(k-cumct[lev-1L], floor(x[lev-1L] / values[lev]))
		for(i in seq(from=maxct, to=0L, by=-1L)){
			ct[lev] <<- i
			x[lev] <<-  x[lev-1L] - i * values[lev]
			cumct[lev] <<- cumct[lev-1L] + i
			if(x[lev] == 0L){
				nsol <<- nsol+1L
				out[, nsol] <<- rep(values, c(k-cumct[lev], ct[2L:lev],rep(0L,nvals-lev)))
			}else if(x[lev] > 0L){
				Recall(lev+1L)
			}else stop('x[levl]<0L')
		}
	}
	#debug(recur.fctn)
	recur.fctn(2L)
	out
}
#restrictedpartitionR(15, 10, 1:7)


restrictedpartitionNR=function(n, k=n, values=seq_len(n))
{  ## non-recursive implementation
	values=c(0L, sort(unique(values[values>0 & values<=n]), decreasing=TRUE))
	nvals=length(values)
	maxv=max(values)
	out=matrix(NA_integer_, k, as.numeric(nparts.atmost.ubound(n,k,maxv)))
	nsol=0L
	ncol.out=ncol(out)

	nextval=c(values[-1L],0L)
	diffvals=values-nextval

	x=rep(n, nvals)
	minct=maxct=ct=rep(0L, nvals)
	ct.allow=rep(k, nvals)

	newLev=TRUE; lev=2L; niter=0L
	repeat{
		niter=niter+1L
		if(newLev){
			while(values[lev]>x[lev-1L] && lev < nvals){
				## fast forward
				ct[lev]=0L
				maxct[lev]=minct[lev]=0L
				ct.allow[lev]=ct.allow[lev-1L]
				x[lev]=x[lev-1L]
				lev=lev+1L
			}
			maxct[lev] = min(ct.allow[lev-1], floor(x[lev-1L] / values[lev]))
			minct[lev] = max(0L, ceiling((x[lev-1L]-ct.allow[lev-1L]*nextval[lev]) / diffvals[lev] ))
			newLev=FALSE
			ct[lev] = maxct[lev] + 1L
		}
		if (ct[lev] > minct[lev]){
			ct[lev] = ct[lev]-1L
		}else{
			lev = lev - 1L
			if(lev == 1L) break else next
		}

		x[lev]  =  x[lev-1L] - ct[lev] * values[lev]
		ct.allow[lev] = ct.allow[lev-1L] - ct[lev]
		if(x[lev] == 0L){
			nsol = nsol+1L
			out[, nsol] = c(rep(values[2:lev], ct[2L:lev]),rep(0L,ct.allow[lev]))
			if(lev==nvals) lev=lev-1L
			if(nsol==ncol.out) break
			#if(nsol%%500==0) cat(niter, nsol, date(), '\n', sep='\t')
		}else if(lev<nvals) {
			lev=lev+1L
			newLev=TRUE
		}

		## fast rewind
		while(ct[lev]==minct[lev] && !newLev && lev>1L) lev=lev-1L
		if(lev==1L) break
	}

	if(nsol<ncol(out)) out=out[, seq_len(nsol), drop=FALSE]
	attr(out, 'niter')=niter
	out
}
#debug(restrictedpartitionNR)
#restrictedpartitionNR(15, 10, 1:7)
#system.time(restrictedpartitionNR(80, 18, 1:9))


restrictedpartitionNRC=function(n, k=n, values=seq_len(n))
{  ## non-recursive implementation
	values=c(0L, sort(unique(as.integer(values[values>0 & values<=n])), decreasing=TRUE))
	nvals=length(values)
	maxv=max(values)
	out=matrix(0L, k, as.numeric(nparts.atmost.ubound(n,k,maxv))) # save time on writing 0's later in C
	nsol=0L
	ncol.out=ncol(out)
	if(ncol.out==0L) return(out)

	nextval=c(values[-1L],0L)
	diffvals=values-nextval

	x=rep(as.integer(n), nvals)
	## do not write the following 3 lines together!
	## they share the same mem location initially
	ct=rep(0L, nvals)
	minct=rep(0L, nvals)
	maxct=rep(0L, nvals)
	ct.allow=rep(as.integer(k), nvals)


	x; ct; minct; maxct; ct.allow;
	values; nextval; diffvals;
	out; ncol.out
	tmpans=.Call(restrParts,
		x, ct, minct, maxct, ct.allow,
		values, nextval, diffvals,
		out, ncol.out
	)
	#browser()

	if(isTRUE(tmpans)) out else tmpans
}
#debug(restrictedpartitionNR)
#restrictedpartitionNR(15, 10, 1:7)
#system.time(restrictedpartitionNR(80, 18, 1:9))


restrictedpartitionZS1=function(n, k=n, upper=n)
{  ## non-recursive implementation
	values=c(0L, as.integer(seq(from=upper, to=1L, by=-1L)))
	nvals=length(values)
	maxv=max(values)
	out=matrix(NA_integer_, k, as.numeric(nparts.atmost.ubound(n,k,maxv)))
	nsol=0L
	ncol.out=ncol(out)
	if(ncol.out==0L) return(out)

	nextval=c(values[-1L],0L)
	diffvals=values-nextval

	x=rep(as.integer(n), nvals)
	## do not write the following 3 lines together!
	## they share the same mem location initially
	ct=rep(0L, nvals)
	minct=rep(0L, nvals)
	maxct=rep(0L, nvals)
	ct.allow=rep(as.integer(k), nvals)


	x; ct; minct; maxct; ct.allow;
	values; nextval; diffvals;
	out; ncol.out
	tmpans=as.vector(.Call(restrParts,
		x, ct, minct, maxct, ct.allow,
		values, nextval, diffvals,
		out, 1L
	))
	startx=out[,1L]
	startx=sort(startx[startx>0], decreasing=TRUE)
	min1st = as.integer(ceiling(n/k))
	startx; min1st
	tmpans2=.Call(ZS1,
		as.integer(n), startx, as.integer(k),
		out, ncol.out, min1st
	)
	#browser()

	if(isTRUE(tmpans2)) out else tmpans2
}
#debug(restrictedpartitionNR)
#restrictedpartitionNR(15, 10, 1:7)
#system.time(restrictedpartitionNR(80, 18, 1:9))

restrictedpartitionRJa=function(n, k=n, lower=0L, upper=n, distinct=FALSE)
{
	num=as.integer(nparts.atmost.ubound(n,k,upper))
	m=as.integer(n);
	l1=as.integer(lower)
	l2=as.integer(upper)
	n=as.integer(k)
	out=matrix(NA_integer_, k, num)
	
	m; l1; l2; n; num; out;
	tmpans= if(distinct) .Call(RJa1 , m, l1, l2, n, num, out) else 
	.Call(RJa0 , m, l1, l2, n, num, out)
	if(isTRUE(tmpans)) out else tmpans
}

restrictedpartitionRJb=function(n, k=n, values=0:n, max.multiplicity = rep(k, length(values)))
{
	num=as.integer(nparts.atmost.ubound(n,k,max(values)))
	r=max(c(length(values), length(max.multiplicity)))
	mu=rep(as.integer(max.multiplicity), length=r)
	m=as.integer(n);
	v=rep(as.integer(values), length=r)
	n=as.integer(k)
	out=matrix(NA_integer_, k, num)
	
	mu; m; values; n; num; out;
	tmpans=.Call(RJb, mu, m, values, n, num, out)
	if(isTRUE(tmpans)) out else tmpans
}

restrictedpartitionRJaSS=function(n, k=n, lower=0L, upper=n, sum.square.upper=n*n)
{
	l2=as.integer(min(n, upper, (n+sqrt((k-1)*(k*sum.square.upper - n*n)))/k))
	ii = floor(n/l2); 
	if(sum.square.upper >= (ii+ii*ii)*l2*l2-2*ii*n*l2+n*n) return( restrictedpartitionRJa(n, k, lower, l2))
	num=as.integer(nparts.atmost.ubound(n,k,l2))
	m=as.integer(n);
	l1=as.integer(lower)
	n=as.integer(k)
	ss=as.integer(sum.square.upper)
	out=matrix(NA_integer_, k, num)
	
	m; l1; l2; n; ss; num; out;
	tmpans=.Call(RJa0ss, m, l1, l2, n, ss, num, out)
	if(isTRUE(tmpans)) out else tmpans
}

partsYKKN=function(n, worstCaseO1=FALSE)
{
	worstCaseO1=as.logical(worstCaseO1)
	#stopifnot(method==1L || method==2L)
	n=as.integer(n);
	out=matrix(0L, n, partitions::P(n));

	n; out; worstCaseO1
	tmpans=.Call(ykknAllParts, n, out, worstCaseO1);
	if(isTRUE(tmpans)) out else browser()
}

restrictedpartsYKKN=function(n, k, worstCaseO1=FALSE)
{
	worstCaseO1=as.logical(worstCaseO1)
	#stopifnot(method==1L || method==2L)
	n=as.integer(n);
	k=min(n, as.integer(k))
	out=matrix(0L, k, as.integer(nparts.atmost(n,k)));

	n; k; out; worstCaseO1
	tmpans=.Call(ykknAtMostKParts, n, k, out, worstCaseO1);
	if(isTRUE(tmpans)) out else browser()
}
