##checking input
test_that("Input check", {

    test1 = matrix(rnorm(100*4), nrow=100, ncol=4) #test input is in matrix form

    #check input has to be GRangesList
    checkException(run.cin.chr(grl.seg = test1, thr.gain=2.25, thr.loss=1.75, V.def=3, V.mode="sum"), silent = TRUE)
})
