library(SBayesRC)
library(data.table)
library(dplyr)
library(tidyverse)
library(logger)

Args <- commandArgs()
block_in <- Args[6]

logger.begin <- function(outPrefix, log2file){
    bLoggerTOKEN <<- log2file
    if(bLoggerTOKEN){
        zzLOGGER <<- file(paste0(outPrefix, ".log"), "wt")
        message("The messages are redirected to ", outPrefix, ".log")
        message("  No output here...")
        sink(zzLOGGER)
        sink(zzLOGGER, type="message")
    }
}

logger.end <- function(){
    if(bLoggerTOKEN){
        sink()
        sink(type="message")
        close(zzLOGGER)
        message("  Output reactivated.")
    }
}


LDstep3 = function(outDir, blockIndex, thresh=0.995, log2file=FALSE){
    message("Step3: perform eigen decomposition on LD matrix")
    message("This step must run after step2")
    output = outDir
    dt = fread(file.path(output, "ldm.info"))

    curLine = as.integer(blockIndex)
    act_block = dt[curLine]$Block
    outFile = file.path(output, paste0("block", act_block, ".eigen.bin"))

    if(file.exists(paste0(outFile))){
        message("The eigen decomposition was done on this block before.")
        message("Block index: ", curLine)
        return
    }

    logger.begin(outFile, log2file)
    if(log2file){
        message("Step3: perform eigen decomposition on LD matrix")
        message("This step must run after step2")
    }

    t_begin = proc.time()

    filename = file.path(output, paste0("b", act_block, ".ldm.full"))
    message(" Block index: ", curLine)
    message(" Block number: ", act_block)
    message(" Output: ", outFile)

    message(" Read LD information ", filename)

    info = fread(paste0(filename, ".info"))
    nMarker = as.numeric(nrow(info))

    ldfile = file(paste0(filename, ".bin"), "rb")
    R = readBin(ldfile, n = nMarker * nMarker, what=numeric(0), size=4)
    dim(R) = c(nMarker, nMarker)
    close(ldfile)

    message(" Start eigen decomposition: m = ", nrow(R))

    eig = eigen(R, symmetric=TRUE)

    rm(R)
    gc()

    lambda = eig$values
    m = length(lambda)

    selected = (lambda > 0)
    if(length(rle(selected)$values) > 2){
        stop("strange eigen pattern")
    }
    wholeLambda = lambda[selected]

    message("  Cut to variance threshold: ", thresh)
    lambdaSum = sum(wholeLambda)
    k = which(cumsum(wholeLambda/lambdaSum) >= thresh)[1]

    message(" k: ", k, ", m: ", m)
    cFile = file(outFile, "wb")
    writeBin(m, cFile, size=4)
    writeBin(k, cFile, size=4)
    writeBin(lambdaSum, cFile, size=4)
    writeBin(thresh, cFile, size=4)
    message(" Write binary...")
    writeBin(wholeLambda[1:k], cFile, size=4)

    for(idx in 1:k){
        writeBin(eig$vectors[, idx], cFile, size=4)
    }

    close(cFile)

    message("Done.")
    print(proc.time() - t_begin)
    logger.end()
    if(log2file){
        print(proc.time() - t_begin)
        message("Done.")
    }
}

LDstep3(outDir="/scratch/gy71651/Huifang/SBayesNew/O3_LD/", blockIndex=block_in, log2file=TRUE)
