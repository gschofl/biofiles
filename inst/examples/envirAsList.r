#---------------------------------------------------------------------------
# ENVIRONMENTS
#---------------------------------------------------------------------------

envir <- new.env()
envir$a <- new.env()
envir$a$a.1 <- new.env()
envir$a$a.1$a.1.1 <- new.env()
envir$a$a.1$a.1.1$a.1.1.1 <- NA
envir$a$a.1$a.1.2 <- 5
envir$a$a.2 <- new.env()
envir$a$a.2$a.2.1 <- list()
envir$b <- NULL

envirAsList(src=envir)

flatten(envir)
flatten(envir, start.after=1)
flatten(envir, start.after=1, .do.debug=TRUE)
flatten(envir, start.after=2)
flatten(envir, start.after=3)
flatten(envir, start.after=4)

flatten(envir, stop.at=1)
flatten(envir, stop.at=2)
flatten(envir, stop.at=3)
flatten(envir, stop.at=4)

flatten(envir, start.after=1, stop.at=1)
flatten(envir, start.after=1, stop.at=2)
flatten(envir, start.after=1, stop.at=3)
flatten(envir, start.after=1, stop.at=4)

# /ENVIRONMENTS ----------
