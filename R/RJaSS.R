if(FALSE){
# this is the R version for debugging purpose only
RJa0ssR=function(m, l1=0L, l2=m, n=m, ss=m*m)
{
	m0=m=as.integer(m);
	l1=as.integer(l1); 
	l2=as.integer(l2);
	n=as.integer(n);
	ss0=ss=as.integer(ss);
	nsols=as.integer(nparts.atmost.ubound(m,n,l2))
	out=matrix(NA_integer_, n, nsols)
	
	l1sq = l1 * l1;  l12= 2L * l1;
	x=rep(l1, n)

	num = 0L; m = m - n * l1 ; l2 = l2 - l1;  ##// diff bounds
	ss = ss - n * l1sq; 
	if (m>=0L && m <= n*l2 && ss>=0L && ss <= n * l2 * l2){
		
		i=1L; ##// l2 -= z*(n-1);
		repeat{
##//a1:
			status=list(x=x,m=m,ss=ss,i=i,l2=l2)
			##// keep allocating current allowed upper bound
			while (m > (toAdd <-as.integer(min(l2, sqrt(l1sq+ss)-l1))) && toAdd > 0L && i<n) 
			{##// toAdd holds allowed upper bound
			##//while (ss >= (l2sq = l2 * (l2 + 2 * l1)) && m > l2)
				m = m - toAdd; ss = ss - toAdd * ( toAdd + l12); 
				x[i] = l1 + toAdd; if(ss!=ss0-sum(x[1:i]^2)) browser()
				i=i+1L
			}
			##//bool found=false;
			if (toAdd != 0L && i<=n && m<=toAdd) {
				x[i] = l1 + m; num =num+1L;  ##// last part;
				##//found = true;
				if(num==18L)browser()
				if(sum(x)!=m0 || sum(x*x) >ss0) browser()
				out[,num]=x
				
				if (num == nsols) break; ##// goto end;
			
				if (i<n && m>1L) { #// possible to increase # of parts?
								  ##// if so, do it;
								  #// but after this, no way to extend it w/o affecting decreasing order -- have to backtrack then.
					##// note that this never decrease ss: 
					##// original sum of squares at location i and i+1 is 
					##//		(l1+toAdd)^2 + l1^2; 
					##// new sum of squares at these 2 locations is 
					##//		(l1+toAdd-1)^2 + (l1+1)^2 = (original ss) -2*(toAdd-1).
					ss = ss - (x[i]-1L)*(x[i]-1L); #// ss=ss+x[i]^2-(x[i]-1)^2
					m = 1L; x[i]=x[i]-1L; 
					if(ss!=ss0-sum(x[1:i]^2)) browser()
					i=i+1L; x[i] = l1 + 1L; 
					num = num + 1L;
					if(num==18L)browser();
					out[,num]=x;
					#for(int tmp=1; tmp<=n; ++tmp) *(out++) = x[tmp];
					if (num == nsols) break; #//goto end;
				}
			}
			j=i-1L
			while(j>=1L){
				if(ss!=ss0-sum(x[1:j]^2))browser()
			#for (j=i-1; j>=1; --j) { // backtracking
				ss = ss + 2L*x[j] - 1L; #// ss + x[j]^2 - (x[j]-1)^2
				if(ss!=ss0-sum(x[seq_len(j-1)]^2)-(x[j]-1)^2)browser()
				l2 = x[j] - l1 - 1L; m=m+1L;  #// go back to location j and consider the next smaller value;
				if (m <= (n-j) * l2 #// at most n-j full terms to allocate after location j
				) { #// is it possible for all >j locations to get a solution?
					#// this is not sufficient condition, but necessary.
					x[j] = l1 + l2;
					break; #// goto a1;
				}else{  #// otherwise, keep going back
					status=list(ss=ss, j=j, x=x, m=m, i=i)
					#ss = ss+ (x[j]-1L)*(x[j]-1L) - l1sq; 
					ss = ss -2L*x[j] + 1L ## recover ss since x[j] did not change
					ss = ss + x[j]*x[j] - l1sq ## ss before loc j if all rest are l1
					m = m+l2; x[i] = l1; i = j; #// i always points to j's next location.  At the end: i=1.
					if(ss!=ss0-sum(x[seq_len(i-1)]^2))browser()
				}
				j=j-1L
			}
			if(i <= 1L) break
		}
	}

	out[,seq_len(num), drop=FALSE]
}
#str(RJa0ssR(75,n=10,ss=1000))
}
