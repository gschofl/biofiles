#-------------------------------------------------------------------------------
# SIMPLE LISTS
#-------------------------------------------------------------------------------

src <- list(DF=data.frame(A=c(1,2)), vec=c("a", "b")) 
src <- list(src,src)

flatten(src) 

# /SIMPLE LISTS ----------

#-------------------------------------------------------------------------------
# NESTED LISTS
#-------------------------------------------------------------------------------

src <- list(a=list(a.1=list(a.1.1=list(a.1.1.1=NA), a.1.2=5), 
                   a.2=list(a.2.1=list())), b=NULL)

flatten(src)
flatten(src, start.after=1)
flatten(src, start.after=1, .do.debug=TRUE)
flatten(src, start.after=2)
flatten(src, start.after=3)
flatten(src, start.after=4)

flatten(src, stop.at=1)
flatten(src, stop.at=2)
flatten(src, stop.at=3)
flatten(src, stop.at=4)

flatten(src, start.after=1, stop.at=1)
flatten(src, start.after=1, stop.at=2)
flatten(src, start.after=1, stop.at=3)
flatten(src, start.after=1, stop.at=4)
flatten(src, start.after=2, stop.at=4)

# /NESTED LISTS ----------
