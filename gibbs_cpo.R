# FUNCTIONS ---------------------------------------------------------------

# Full conditionals -------------------------------------------------------

full_U <- function(i, Vs, Zs, invSigt, invSds){

  mu <- chol2inv(chol(invSigt+Vs%*%invSds%*%t(Vs)))%*%(invSigt%*%mU+Vs%*%invSds%*%Zs[i,])
  sig <- chol2inv(chol(invSigt+Vs%*%invSds%*%t(Vs)))
  rmnorm(1, mu, sig)

}

full_R <- function(i, Ws, Ys, invSig, invSdt){

  phi3 <- 1/sqrt(nrow(Ws))
  mu <- chol2inv(chol(1/phi3*invSig+Ws%*%invSdt%*%t(Ws)))%*%(1/phi3*invSig%*%mR+Ws%*%invSdt%*%Ys[i,])
  sig <- chol2inv(chol(1/phi3*invSig+Ws%*%invSdt%*%t(Ws)))
  rmnorm(1, mu, sig)

}

full_V <- function(j, Us, Zs, invUpst, invSns){

  mu <- chol2inv(chol(invUpst+t(Us)%*%invSns%*%Us))%*%(invUpst%*%mV+t(Us)%*%invSns%*%Zs[,j])
  sig <- chol2inv(chol(invUpst+t(Us)%*%invSns%*%Us))
  rmnorm(1, mu, sig)

}

full_W <- function(j, Rs, Ys, invUps, invSnt){

  omega3 <- 1/sqrt(ncol(Rs))
  mu <- chol2inv(chol(1/omega3*invUps+t(Rs)%*%invSnt%*%Rs))%*%(1/omega3*invUps%*%mW+t(Rs)%*%invSnt%*%Ys[,j])
  sig <- chol2inv(chol(1/omega3*invUps+t(Rs)%*%invSnt%*%Rs))
  rmnorm(1, mu, sig)

}

# Gibbs function ----------------------------------------------------------

gibbs_cpo <- function(nIter, Us, Rs, Vs, Ws, Zs, Ys, sigma2, tau2, type_data="binary"){

  # Zs=Z; Ys=Y; Us=U; Ws=W; Rs=R; Vs=V
  set.seed(125)
  n <- nrow(Zs)
  p <- ncol(Zs)
  d <- ncol(Us)

  c_row = 1:n
  c_col = 1:p

  l_v <- c(max(c_col))
  k_v <- c(max(c_row))

  phi1 <- phi3 <- 1/sqrt(d)
  omega1 <- omega3 <- 1/sqrt(d)

  invSds <- 1/sigma2*diag(1,p)
  invSns <- 1/sigma2*diag(1,n)
  invSdt <- 1/tau2*diag(1,p)
  invSnt <- 1/tau2*diag(1,n)
  Sigma <- diag(1, d)
  Sigma_tilde <- phi1*Sigma
  Upsilon <- diag(1, d)
  Upsilon_tilde <- omega1*Upsilon
  invSigt <- chol2inv(chol(Sigma_tilde))
  invSig <- chol2inv(chol(Sigma))
  invUpst <- chol2inv(chol(Upsilon_tilde))
  invUps <- chol2inv(chol(Upsilon))

  temp_cpo <- matrix(0, nrow=n, ncol=p)

  not_def <- c(Inf, -Inf, NaN, NA)

  partition_list <- vector(mode = "list", length = nIter)

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nIter,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar


  for(r in 1:nIter){

    pb$tick()

    for(i in 1:n){

      cons_i <- -1/2*log(det(sigma2*tau2*diag(1,p)))-1/2*log(phi3)-
        1/2*log(det(Sigma%*%Sigma_tilde))-1/2*log(det((invSigt+Vs%*%invSds%*%t(Vs))%*%(1/phi3*invSig+Ws%*%invSdt%*%t(Ws))))-
        1/2*sum(diag(invSds%*%Zs[i,]%*%t(Zs[i,])+invSdt%*%Ys[i,]%*%t(Ys[i,])))-
        1/2*sum(diag(1/phi3*invSig%*%mR%*%t(mR)+invSigt%*%mU%*%t(mU)))-
        1/2*sum(diag((Ws%*%invSdt%*%Ys[i,]+1/phi3*invSig%*%mR)%*%t((Ws%*%invSdt%*%Ys[i,]+1/phi3*invSig%*%mR))%*%chol2inv(chol(Ws%*%invSdt%*%t(Ws)+1/phi3*invSig))))-
        1/2*sum(diag((Vs%*%invSds%*%Zs[i,]+invSigt%*%mU)%*%t((Vs%*%invSds%*%Zs[i,]+invSigt%*%mU))%*%chol2inv(chol(Vs%*%invSds%*%t(Vs)+invSigt))))

      log_pnew_r <- log(Mh) + cons_i

      c_row_i <- c_row[-i]
      new_r <- rep(0, length(c_row_i))
      log_p_r <- rep(NA, length(unique(c_row_i)))

      for(k in unique(c_row_i)){

        nk <- sum(c_row_i==k)

        if(new_r[which(c_row_i==k)[1]]==0){

          if(nk>1){
            Uk <- Us[which(c_row_i==k),][1,]
            Rk <- Rs[which(c_row_i==k),][1,]
          }
          if(nk==1){
            Uk <- Us[which(c_row_i==k),]
            Rk <- Rs[which(c_row_i==k),]
          }

          log_p_r[k] <- log(nk)-
            1/2*sum(diag(invSds%*%as.matrix(Zs[i,]- t(Vs)%*%Uk)%*%as.matrix(t(Zs[i,]- t(Vs)%*%Uk))))-
            1/2*sum(diag(invSdt%*%as.matrix(Ys[i,]- t(Ws)%*%Rk)%*%as.matrix(t(Ys[i,]- t(Ws)%*%Rk))))

          new_r[which(c_row_i==k)] <- 1

        }

      }

      log_p_r <- log_p_r[!is.na(log_p_r)]

      if(log_pnew_r%in%not_def){log_pnew_r <- min(log_p_r)}

      maxr <- max(c(log_pnew_r, log_p_r))

      pnew_r <- exp(log_pnew_r - maxr)

      p_r <- exp(log_p_r - maxr)

      q_r <- sapply(c(pnew_r,p_r), function(x) x/sum(c(p_r,pnew_r))) # normalized probs

      pos <- match(sample(q_r,size=1, prob=q_r), q_r)

      if(pos!=1){

        npos <- sum(c_row==(pos-1))

        if(npos>1){

          Us[i,] <- Us[which(c_row==(pos-1)),][1,]
          Rs[i,] <- Rs[which(c_row==(pos-1)),][1,]


        }
        if(npos==1){
          Us[i,] <- Us[which(c_row==(pos-1)),]
          Rs[i,] <- Rs[which(c_row==(pos-1)),]

        }

        c_row[i] <- (pos-1)

        for(num in 1:max(c_row)){

          if(!(num %in% c_row)){

            c_row[which(c_row==(num+1))] <- c_row[which(c_row==(num+1))]-1
          }

        }

      }
      if(pos==1){

        n_i <- sum(c_row==c_row[i])

        if(n_i!=1){

          c_row[i] <- max(c_row)+1

        }
        Us[i,] <- full_U(i, Vs, Zs, invSigt, invSds)
        Rs[i,] <- full_R(i, Ws, Ys, invSig, invSdt)
      }

    }
    for(num in 1:max(c_row)){

      new_max <- max(c_row)

      if(num<=new_max){

        if(!(num %in% c_row)){

          el <- min(c_row[which(c_row>num)])

          pos <- which(c_row>num)

          c_row[pos] <- c_row[pos]-(el-num)

        }
      }

    }

    # acceleration step

    for(k in c_row){

      nk <- sum(c_row==k)

      new_r <- rep(0, length(c_row))

      if(new_r[k]==0){

        Zk <- colMeans(matrix(Zs[which(c_row==k),], ncol=p))

        Yk <- colMeans(matrix(Ys[which(c_row==k),], ncol=p))

        Us[which(c_row==k),] <- rmnorm(1, chol2inv(chol(invSigt+nk*Vs%*%invSds%*%t(Vs)))%*%(invSigt%*%mU+nk*Vs%*%invSds%*%Zk), chol2inv(chol(invSigt+nk*Vs%*%invSds%*%t(Vs))))
        Rs[which(c_row==k),] <- rmnorm(1, chol2inv(chol(1/phi3*invSig+nk*Ws%*%invSdt%*%t(Ws)))%*%(1/phi3*invSig%*%mR+nk*Ws%*%invSdt%*%Yk), chol2inv(chol(1/phi3*invSig+nk*Ws%*%invSdt%*%t(Ws))))

        new_r[k] <- 1

      }

      if(nk>1){
        Uk <- Us[which(c_row==k),][1,]
        Rk <- Rs[which(c_row==k),][1,]
      }
      if(nk==1){
        Uk <- Us[which(c_row==k),]
        Rk <- Rs[which(c_row==k),]
      }
    }

    z <- sqrt(1/(sqrt(d)*var(as.vector(Us))))
    Us <- z*Us
    Vs <- Vs/z
    M <- varimax(Us)$rotmat
    Us <- Us%*%M
    Vs <- t(M)%*%Vs

    z <- sqrt(1/(sqrt(d)*var(as.vector(Rs))))
    Rs <- z*Rs
    Ws <- Ws/z
    M <- varimax(Rs)$rotmat
    Rs <- Rs%*%M
    Ws <- t(M)%*%Ws

    for(j in 1:p){

      cons_j <- -1/2*log(det(sigma2*tau2*diag(1,n)%*%diag(1,n)))-1/2*log(omega3)-
        1/2*log(det(Upsilon%*%Upsilon_tilde))-1/2*log(det((invUpst+t(Us)%*%invSns%*%Us)%*%(1/omega3*invUps+t(Rs)%*%invSnt%*%Rs)))-
        1/2*sum(diag(invSns%*%Zs[,j]%*%t(Zs[,j])+invSnt%*%Ys[,j]%*%t(Ys[,j])))-
        1/2*sum(diag(1/omega3*invUps%*%mW%*%t(mW)+invUpst%*%mV%*%t(mV)))-
        1/2*sum(diag((t(Rs)%*%invSnt%*%Ys[,j]+1/omega3*invUps%*%mW)%*%t((t(Rs)%*%invSnt%*%Ys[,j]+1/omega3*invUps%*%mW))%*%chol2inv(chol(t(Rs)%*%invSnt%*%Rs+1/omega3*invUps))))-
        1/2*sum(diag((t(Us)%*%invSns%*%Zs[,j]+invUpst%*%mV)%*%t((t(Us)%*%invSns%*%Zs[,j]+invUpst%*%mV))%*%chol2inv(chol(t(Us)%*%invSns%*%Us+invUpst))))

      log_pnew_c <- log(Mf) + cons_j

      c_col_j <- c_col[-j]
      new_c <- rep(0, length(c_col_j))
      log_p_c <- rep(NA, length(unique(c_col_j)))

      for(l in unique(c_col_j)){

        nl <- sum(c_col_j==l)

        if(new_c[which(c_col_j==l)[1]]==0){

          if(nl>1){
            Vl <- Vs[,which(c_col_j==l)][,1]
            Wl <- Ws[,which(c_col_j==l)][,1]
          }else{
            Vl <- Vs[,which(c_col_j==l)]
            Wl <- Ws[,which(c_col_j==l)]
          }

          log_p_c[l] <- log(nl)-
            1/2*sum(diag(invSns%*%as.matrix(Zs[,j]- Us%*%Vl)%*%as.matrix(t(Zs[,j]- Us%*%Vl))))-
            1/2**sum(diag(invSnt%*%as.matrix(Ys[,j]- Rs%*%Wl)%*%as.matrix(t(Ys[,j]- Rs%*%Wl))))

          new_c[which(c_col_j==l)] <- 1
        }
      }

      log_p_c <- log_p_c[!is.na(log_p_c)]

      if(log_pnew_c%in%not_def){log_pnew_c <- min(log_p_c)}

      maxc <- max(c(log_pnew_c, log_p_c))

      pnew_c <- exp(log_pnew_c - maxc)

      p_c <- exp(log_p_c - maxc)

      q_c <- sapply(c(pnew_c,p_c), function(x) x/sum(c(p_c,pnew_c)))

      pos <- match(sample(q_c,size=1, prob=q_c), q_c)

      if(pos!=1){

        if(sum(c_col==(pos-1))>1){

          Vs[,j] <- Vs[,which(c_col==(pos-1))][,1]
          Ws[,j] <- Ws[,which(c_col==(pos-1))][,1]

        }else{

          Vs[,j] <- Vs[,which(c_col==(pos-1))]
          Ws[,j] <- Ws[,which(c_col==(pos-1))]
        }

        c_col[j] <- (pos-1)

        for(num in 1:max(c_col)){

          if(!(num %in% c_col)){

            c_col[which(c_col==(num+1))] <- c_col[which(c_col==(num+1))]-1
          }

        }

      }else{ # sample a new value

        n_j <- sum(c_col==c_col[j])

        if(n_j!=1){

          c_col[j] <- max(c_col)+1

        }

        Vs[,j] <- full_V(j, Us, Zs, invUpst, invSns)
        Ws[,j] <- full_W(j, Rs, Ys, invUps, invSnt)
      }
    }

    for(l in c_col){

      nl <- sum(c_col==l)

      new_c <- rep(0, length(c_col))

      if(new_c[l]==0){

        Zl <- rowMeans(matrix(Zs[,which(c_col==l)], nrow=n))

        Yl <- rowMeans(as.matrix(Ys[,which(c_col==l)], nrow=n))

        Vs[,which(c_col==l)] <- rmnorm(1, chol2inv(chol(invUpst+t(Us)%*%invSns%*%Us))%*%(invUpst%*%mV+nl*t(Us)%*%invSns%*%Zl), chol2inv(chol(invUpst+nl*t(Us)%*%invSns%*%Us)))
        Ws[,which(c_col==l)] <- rmnorm(1, chol2inv(chol(1/omega3*invUps+t(Rs)%*%invSnt%*%Rs))%*%(1/omega3*invUps%*%mW+nl*t(Rs)%*%invSnt%*%Yl), chol2inv(chol(1/omega3*invUps+nl*t(Rs)%*%invSnt%*%Rs)))


        new_c[l] <- 1
      }

      if(nl>1){
        Vl <- Vs[,which(c_col==l)][,1]
        Wl <- Ws[,which(c_col==l)][,1]
      }else{
        Vl <- Vs[,which(c_col==l)]
        Wl <- Ws[,which(c_col==l)]
      }


    }

    for(num in 1:max(c_col)){

      new_max <- max(c_col)

      if(num<=new_max){

        if(!(num %in% c_col)){

          el <- min(c_row[which(c_col>num)])

          pos <- which(c_col>num)

          c_col[pos] <- c_col[pos]-(el-num)

        }
      }

    }

    m <- sqrt(1/(sqrt(d)*var(as.vector(Vs))))
    Vs <- Vs*m
    Us <-(Us/m)
    M <- varimax(t(Vs))$rotmat
    Vs <- t(M)%*%(Vs)
    Us <- Us%*%M

    m <- sqrt(1/(sqrt(d)*var(as.vector(Ws))))
    Ws <- Ws*m
    Rs <-(Rs/m)
    M <- varimax(t(Ws))$rotmat
    Ws <- t(M)%*%(Ws)
    Rs <- Rs%*%M

    for(i in 1:n){
      for(j in 1:p){

        if(type_data=="ordinal"){if(is.na(data[i,j])==F){

          if(data[i,j]==1){

            Zs[i,j] <- rtnorm(1, lower=-Inf, upper=-0.5, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))




          }
          if(data[i,j]==2){


            Zs[i,j] <- rtnorm(1, lower =-0.5, upper =0.5, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))



          }
          if(data[i,j]==3){


            Zs[i,j] <- rtnorm(1, lower=0.5, upper=Inf, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))




          }

          Ys[i,j] <- rtnorm(1, lower=0, upper=Inf, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }else{

          Zs[i,j] <- rnorm(1, Us[i,]%*%Vs[,j], sqrt(sigma2))
          Ys[i,j] <- rtnorm(1, lower=-Inf, upper=0, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }
        }

        if(type_data=="binary"){if(is.na(data[i,j])==F){

          if(data[i,j]==1){
            Zs[i,j] <- rtnorm(1, lower=0, upper=Inf, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))

          }
          if(data[i,j]==0){

            Zs[i,j] <- rtnorm(1, lower =-Inf, upper =0, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))


          }

          Ys[i,j] <- rtnorm(1, lower=0, upper=Inf, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }else{

          Zs[i,j] <- rnorm(1, Us[i,]%*%Vs[,j], sqrt(sigma2))
          Ys[i,j] <- rtnorm(1, lower=-Inf, upper=0, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }}

      }


    }

    Zs <- scale(Zs)
    Ys <- scale(Ys)

    if(type_data=="ordinal"){
      for(i in 1:n){
        for(j in 1:p){

          if(is.na(data[i,j])==F){

            if(data[i,j]==1){
              newz <- -log(pnorm(-0.5, Us[i,]%*%Vs[,j], sqrt(sigma2))+1)
            }
            if(data[i,j]==2){
              newz <- -log(pnorm(0.5, Us[i,]%*%Vs[,j], sqrt(sigma2))-pnorm(-0.5, Us[i,]%*%Vs[,j], sqrt(sigma2))+1)

            }
            if(data[i,j]==3){
              newz <- -log(1-pnorm(0.5, Us[i,]%*%Vs[,j], sqrt(sigma2))+1)

            }
            newy <- -log(1-pnorm(0, Rs[i,]%*%Ws[,j], sqrt(tau2))+1)

          }else{
            newy <- -log(pnorm(0, Rs[i,]%*%Ws[,j], sqrt(tau2))+1)
            newz <- -log(1+1)

          }

          temp_cpo[i,j] <- temp_cpo[i,j]+exp(newz+newy)

        }
      }
    }
    if(type_data=="binary"){

      for(i in 1:n){
        for(j in 1:p){

          if(is.na(data[i,j])==F){

            if(data[i,j]==1){
              newz <- -log(1-pnorm(0, Us[i,]%*%Vs[,j], sqrt(sigma2))+1)
            }
            if(data[i,j]==0){
              newz <- -log(pnorm(0, Us[i,]%*%Vs[,j], sqrt(sigma2))+1)

            }
            newy <- -log(1-pnorm(0, Rs[i,]%*%Ws[,j], sqrt(tau2))+1)

          }else{
            newy <- -log(pnorm(0, Rs[i,]%*%Ws[,j], sqrt(tau2))+1)
            newz <- 0
          }

          temp_cpo[i,j] <- temp_cpo[i,j]+exp(newz+newy)
        }
      }

    }

    k_v <- c(k_v, max(c_row))
    l_v <- c(l_v, max(c_col))

    lpml <- n*p*log(nIter) - sum(rowSums(log(temp_cpo)))

    partition_list[r] <- list(list(c_row, c_col))

  }


  list_data <- list(c_row, c_col,
                    partition_list,
                    k_v, l_v, lpml

  )
  names(list_data) <- c("c_row",
                        "c_col",
                        "partition_list",
                        "k_v", "l_v", "lpml")
  return(list_data)
}

gibbs <- function(nIter, Us, Rs, Vs, Ws, Zs, Ys, sigma2, tau2, type_data="binary"){

  set.seed(125)
  n <- nrow(Zs)
  p <- ncol(Zs)
  d <- ncol(Us)

  c_row = 1:n
  c_col = 1:p

  l_v <- c(max(c_col))
  k_v <- c(max(c_row))

  phi1 <- phi3 <- 1/sqrt(d)
  omega1 <- omega3 <- 1/sqrt(d)

  invSds <- 1/sigma2*diag(1,p)
  invSns <- 1/sigma2*diag(1,n)
  invSdt <- 1/tau2*diag(1,p)
  invSnt <- 1/tau2*diag(1,n)
  Sigma <- diag(1, d)
  Sigma_tilde <- phi1*Sigma
  Upsilon <- diag(1, d)
  Upsilon_tilde <- omega1*Upsilon
  invSigt <- chol2inv(chol(Sigma_tilde))
  invSig <- chol2inv(chol(Sigma))
  invUpst <- chol2inv(chol(Upsilon_tilde))
  invUps <- chol2inv(chol(Upsilon))

  not_def <- c(Inf, -Inf, NaN, NA)

  partition_list <- vector(mode = "list", length = nIter)

  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nIter,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar


  for(r in 1:nIter){

    pb$tick()

    for(i in 1:n){

      cons_i <- -1/2*log(det(sigma2*tau2*diag(1,p)))-1/2*log(phi3)-
        1/2*log(det(Sigma%*%Sigma_tilde))-1/2*log(det((invSigt+Vs%*%invSds%*%t(Vs))%*%(1/phi3*invSig+Ws%*%invSdt%*%t(Ws))))-
        1/2*sum(diag(invSds%*%Zs[i,]%*%t(Zs[i,])+invSdt%*%Ys[i,]%*%t(Ys[i,])))-
        1/2*sum(diag(1/phi3*invSig%*%mR%*%t(mR)+invSigt%*%mU%*%t(mU)))-
        1/2*sum(diag((Ws%*%invSdt%*%Ys[i,]+1/phi3*invSig%*%mR)%*%t((Ws%*%invSdt%*%Ys[i,]+1/phi3*invSig%*%mR))%*%chol2inv(chol(Ws%*%invSdt%*%t(Ws)+1/phi3*invSig))))-
        1/2*sum(diag((Vs%*%invSds%*%Zs[i,]+invSigt%*%mU)%*%t((Vs%*%invSds%*%Zs[i,]+invSigt%*%mU))%*%chol2inv(chol(Vs%*%invSds%*%t(Vs)+invSigt))))

      log_pnew_r <- log(Mh) + cons_i

      c_row_i <- c_row[-i]
      new_r <- rep(0, length(c_row_i))
      log_p_r <- rep(NA, length(unique(c_row_i)))

      for(k in unique(c_row_i)){

        nk <- sum(c_row_i==k)

        if(new_r[which(c_row_i==k)[1]]==0){

          if(nk>1){
            Uk <- Us[which(c_row_i==k),][1,]
            Rk <- Rs[which(c_row_i==k),][1,]
          }
          if(nk==1){
            Uk <- Us[which(c_row_i==k),]
            Rk <- Rs[which(c_row_i==k),]
          }

          log_p_r[k] <- log(nk)-
            1/2*sum(diag(invSds%*%as.matrix(Zs[i,]- t(Vs)%*%Uk)%*%as.matrix(t(Zs[i,]- t(Vs)%*%Uk))))-
            1/2*sum(diag(invSdt%*%as.matrix(Ys[i,]- t(Ws)%*%Rk)%*%as.matrix(t(Ys[i,]- t(Ws)%*%Rk))))

          new_r[which(c_row_i==k)] <- 1

        }

      }

      log_p_r <- log_p_r[!is.na(log_p_r)]

      if(log_pnew_r%in%not_def){log_pnew_r <- min(log_p_r)}

      maxr <- max(c(log_pnew_r, log_p_r))

      pnew_r <- exp(log_pnew_r - maxr)

      p_r <- exp(log_p_r - maxr)

      q_r <- sapply(c(pnew_r,p_r), function(x) x/sum(c(p_r,pnew_r))) # normalized probs

      pos <- match(sample(q_r,size=1, prob=q_r), q_r)

      if(pos!=1){

        npos <- sum(c_row==(pos-1))

        if(npos>1){

          Us[i,] <- Us[which(c_row==(pos-1)),][1,]
          Rs[i,] <- Rs[which(c_row==(pos-1)),][1,]


        }
        if(npos==1){
          Us[i,] <- Us[which(c_row==(pos-1)),]
          Rs[i,] <- Rs[which(c_row==(pos-1)),]

        }

        c_row[i] <- (pos-1)

        for(num in 1:max(c_row)){

          if(!(num %in% c_row)){

            c_row[which(c_row==(num+1))] <- c_row[which(c_row==(num+1))]-1
          }

        }

      }
      if(pos==1){

        n_i <- sum(c_row==c_row[i])

        if(n_i!=1){

          c_row[i] <- max(c_row)+1

        }
        Us[i,] <- full_U(i, Vs, Zs, invSigt, invSds)
        Rs[i,] <- full_R(i, Ws, Ys, invSig, invSdt)
      }

    }
    for(num in 1:max(c_row)){

      new_max <- max(c_row)

      if(num<=new_max){

        if(!(num %in% c_row)){

          el <- min(c_row[which(c_row>num)])

          pos <- which(c_row>num)

          c_row[pos] <- c_row[pos]-(el-num)

        }
      }

    }

    # acceleration step

    for(k in c_row){

      nk <- sum(c_row==k)

      new_r <- rep(0, length(c_row))

      if(new_r[k]==0){

        Zk <- colMeans(matrix(Zs[which(c_row==k),], ncol=p))

        Yk <- colMeans(matrix(Ys[which(c_row==k),], ncol=p))

        Us[which(c_row==k),] <- rmnorm(1, chol2inv(chol(invSigt+nk*Vs%*%invSds%*%t(Vs)))%*%(invSigt%*%mU+nk*Vs%*%invSds%*%Zk), chol2inv(chol(invSigt+nk*Vs%*%invSds%*%t(Vs))))
        Rs[which(c_row==k),] <- rmnorm(1, chol2inv(chol(1/phi3*invSig+nk*Ws%*%invSdt%*%t(Ws)))%*%(1/phi3*invSig%*%mR+nk*Ws%*%invSdt%*%Yk), chol2inv(chol(1/phi3*invSig+nk*Ws%*%invSdt%*%t(Ws))))

        new_r[k] <- 1

      }

      if(nk>1){
        Uk <- Us[which(c_row==k),][1,]
        Rk <- Rs[which(c_row==k),][1,]
      }
      if(nk==1){
        Uk <- Us[which(c_row==k),]
        Rk <- Rs[which(c_row==k),]
      }
    }

    z <- sqrt(1/(sqrt(d)*var(as.vector(Us))))
    Us <- z*Us
    Vs <- Vs/z
    M <- varimax(Us)$rotmat
    Us <- Us%*%M
    Vs <- t(M)%*%Vs

    z <- sqrt(1/(sqrt(d)*var(as.vector(Rs))))
    Rs <- z*Rs
    Ws <- Ws/z
    M <- varimax(Rs)$rotmat
    Rs <- Rs%*%M
    Ws <- t(M)%*%Ws

    for(j in 1:p){

      cons_j <- -1/2*log(det(sigma2*tau2*diag(1,n)%*%diag(1,n)))-1/2*log(omega3)-
        1/2*log(det(Upsilon%*%Upsilon_tilde))-1/2*log(det((invUpst+t(Us)%*%invSns%*%Us)%*%(1/omega3*invUps+t(Rs)%*%invSnt%*%Rs)))-
        1/2*sum(diag(invSns%*%Zs[,j]%*%t(Zs[,j])+invSnt%*%Ys[,j]%*%t(Ys[,j])))-
        1/2*sum(diag(1/omega3*invUps%*%mW%*%t(mW)+invUpst%*%mV%*%t(mV)))-
        1/2*sum(diag((t(Rs)%*%invSnt%*%Ys[,j]+1/omega3*invUps%*%mW)%*%t((t(Rs)%*%invSnt%*%Ys[,j]+1/omega3*invUps%*%mW))%*%chol2inv(chol(t(Rs)%*%invSnt%*%Rs+1/omega3*invUps))))-
        1/2*sum(diag((t(Us)%*%invSns%*%Zs[,j]+invUpst%*%mV)%*%t((t(Us)%*%invSns%*%Zs[,j]+invUpst%*%mV))%*%chol2inv(chol(t(Us)%*%invSns%*%Us+invUpst))))

      log_pnew_c <- log(Mf) + cons_j

      c_col_j <- c_col[-j]
      new_c <- rep(0, length(c_col_j))
      log_p_c <- rep(NA, length(unique(c_col_j)))

      for(l in unique(c_col_j)){

        nl <- sum(c_col_j==l)

        if(new_c[which(c_col_j==l)[1]]==0){

          if(nl>1){
            Vl <- Vs[,which(c_col_j==l)][,1]
            Wl <- Ws[,which(c_col_j==l)][,1]
          }else{
            Vl <- Vs[,which(c_col_j==l)]
            Wl <- Ws[,which(c_col_j==l)]
          }

          log_p_c[l] <- log(nl)-
            1/2*sum(diag(invSns%*%as.matrix(Zs[,j]- Us%*%Vl)%*%as.matrix(t(Zs[,j]- Us%*%Vl))))-
            1/2**sum(diag(invSnt%*%as.matrix(Ys[,j]- Rs%*%Wl)%*%as.matrix(t(Ys[,j]- Rs%*%Wl))))

          new_c[which(c_col_j==l)] <- 1
        }
      }

      log_p_c <- log_p_c[!is.na(log_p_c)]

      if(log_pnew_c%in%not_def){log_pnew_c <- min(log_p_c)}

      maxc <- max(c(log_pnew_c, log_p_c))

      pnew_c <- exp(log_pnew_c - maxc)

      p_c <- exp(log_p_c - maxc)

      q_c <- sapply(c(pnew_c,p_c), function(x) x/sum(c(p_c,pnew_c)))

      pos <- match(sample(q_c,size=1, prob=q_c), q_c)

      if(pos!=1){

        if(sum(c_col==(pos-1))>1){

          Vs[,j] <- Vs[,which(c_col==(pos-1))][,1]
          Ws[,j] <- Ws[,which(c_col==(pos-1))][,1]

        }else{

          Vs[,j] <- Vs[,which(c_col==(pos-1))]
          Ws[,j] <- Ws[,which(c_col==(pos-1))]
        }

        c_col[j] <- (pos-1)

        for(num in 1:max(c_col)){

          if(!(num %in% c_col)){

            c_col[which(c_col==(num+1))] <- c_col[which(c_col==(num+1))]-1
          }

        }

      }else{ # sample a new value

        n_j <- sum(c_col==c_col[j])

        if(n_j!=1){

          c_col[j] <- max(c_col)+1

        }

        Vs[,j] <- full_V(j, Us, Zs, invUpst, invSns)
        Ws[,j] <- full_W(j, Rs, Ys, invUps, invSnt)
      }
    }

    for(l in c_col){

      nl <- sum(c_col==l)

      new_c <- rep(0, length(c_col))

      if(new_c[l]==0){

        Zl <- rowMeans(matrix(Zs[,which(c_col==l)], nrow=n))

        Yl <- rowMeans(as.matrix(Ys[,which(c_col==l)], nrow=n))

        Vs[,which(c_col==l)] <- rmnorm(1, chol2inv(chol(invUpst+t(Us)%*%invSns%*%Us))%*%(invUpst%*%mV+nl*t(Us)%*%invSns%*%Zl), chol2inv(chol(invUpst+nl*t(Us)%*%invSns%*%Us)))
        Ws[,which(c_col==l)] <- rmnorm(1, chol2inv(chol(1/omega3*invUps+t(Rs)%*%invSnt%*%Rs))%*%(1/omega3*invUps%*%mW+nl*t(Rs)%*%invSnt%*%Yl), chol2inv(chol(1/omega3*invUps+nl*t(Rs)%*%invSnt%*%Rs)))


        new_c[l] <- 1
      }

      if(nl>1){
        Vl <- Vs[,which(c_col==l)][,1]
        Wl <- Ws[,which(c_col==l)][,1]
      }else{
        Vl <- Vs[,which(c_col==l)]
        Wl <- Ws[,which(c_col==l)]
      }


    }

    for(num in 1:max(c_col)){

      new_max <- max(c_col)

      if(num<=new_max){

        if(!(num %in% c_col)){

          el <- min(c_row[which(c_col>num)])

          pos <- which(c_col>num)

          c_col[pos] <- c_col[pos]-(el-num)

        }
      }

    }

    m <- sqrt(1/(sqrt(d)*var(as.vector(Vs))))
    Vs <- Vs*m
    Us <-(Us/m)
    M <- varimax(t(Vs))$rotmat
    Vs <- t(M)%*%(Vs)
    Us <- Us%*%M

    m <- sqrt(1/(sqrt(d)*var(as.vector(Ws))))
    Ws <- Ws*m
    Rs <-(Rs/m)
    M <- varimax(t(Ws))$rotmat
    Ws <- t(M)%*%(Ws)
    Rs <- Rs%*%M

    for(i in 1:n){
      for(j in 1:p){

        if(type_data=="ordinal"){if(is.na(data[i,j])==F){

          if(data[i,j]==1){

            Zs[i,j] <- rtnorm(1, lower=-Inf, upper=-0.5, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))




          }
          if(data[i,j]==2){


            Zs[i,j] <- rtnorm(1, lower =-0.5, upper =0.5, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))



          }
          if(data[i,j]==3){

            Zs[i,j] <- rtnorm(1, lower=0.5, upper=Inf, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))




          }

          Ys[i,j] <- rtnorm(1, lower=0, upper=Inf, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }else{

          Zs[i,j] <- rnorm(1, Us[i,]%*%Vs[,j], sqrt(sigma2))
          Ys[i,j] <- rtnorm(1, lower=-Inf, upper=0, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }
        }

        if(type_data=="binary"){if(is.na(data[i,j])==F){

          if(data[i,j]==1){
            Zs[i,j] <- rtnorm(1, lower=0, upper=Inf, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))

          }
          if(data[i,j]==0){

            Zs[i,j] <- rtnorm(1, lower =-Inf, upper =0, mean = Us[i,]%*%Vs[,j], sqrt(sigma2))


          }

          Ys[i,j] <- rtnorm(1, lower=0, upper=Inf, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }else{

          Zs[i,j] <- rnorm(1, Us[i,]%*%Vs[,j], sqrt(sigma2))
          Ys[i,j] <- rtnorm(1, lower=-Inf, upper=0, Rs[i,]%*%Ws[,j], sqrt(tau2))

        }}

      }


    }

    Zs <- scale(Zs)
    Ys <- scale(Ys)

    Zs <- scale(Zs)
    Ys <- scale(Ys)

    k_v <- c(k_v, max(c_row))
    l_v <- c(l_v, max(c_col))


    partition_list[r] <- list(list(c_row, c_col))

  }

  list_data <- list(c_row, c_col,
                    partition_list,
                    k_v, l_v

  )
  names(list_data) <- c("c_row",
                        "c_col",
                        "partition_list",
                        "k_v", "l_v")
  return(list_data)
}

# Similarity --------------------------------------------------------------

similarity_funct <- function(N, niter, burn = 0, partition_list, dim = 1){

  nIter <- niter
  niter <- niter-burn

  PARTITION=matrix(rep(0,N*niter), nrow=niter)
  for (j in 1:niter){
    index <- j + burn
    PARTITION[j,]=unlist(partition_list[[index]][dim])
  }

  # OPTIMAL PARTITION
  psm<-comp.psm(PARTITION)
  data.VI=minVI(psm,PARTITION,method=("all"),include.greedy=FALSE)
  optimal_clustering<-data.VI$cl[3,]

  # SIMILARITY PLOT
  #
  # reorder_mat <- function(mat) {
  #   dd <- as.dist((1 - mat) / 2)
  #   hc <- hclust(dd)
  #   mat <- mat[hc$order, hc$order]
  # }
  #
  # # Get upper triangle of the correlation matrix
  # get_upper_tri <- function(mat) {
  #   mat[lower.tri(mat)] <- NA
  #   return(mat)
  # }
  # psm <- reorder_mat(psm)
  # upper_tri <- get_upper_tri(psm)
  # melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # gg <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
  #   geom_tile(color = "white") +
  #   scale_fill_viridis(option = "viridis", name = "Posterior\nsimilarity", direction=-1, na.value = "white", limits = c(0.00 ,1.00)) +
  #   theme_minimal() +
  #   coord_fixed() +
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.title.y = element_blank(),
  #     panel.grid.major = element_blank(),
  #     panel.border = element_blank(),
  #     panel.background = element_blank(),
  #     axis.ticks = element_blank(),
  #     legend.justification = c(1, 0),
  #     legend.position = c(0.6, 0.7),
  #     legend.direction = "horizontal"
  #   ) +
  #   guides(fill = guide_colorbar(
  #     barwidth = 7,
  #     barheight = 1,
  #     title.position = "top",
  #     title.hjust = 0.5
  #   ))
  #
  # print(gg)
  return(optimal_clustering)

}

