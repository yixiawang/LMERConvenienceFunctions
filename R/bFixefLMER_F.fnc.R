bfFixefLMER_F.fnc<-function(model, item = FALSE, alpha = 0.05, llrt = FALSE, 
          prune.ranefs = TRUE, p.value = "upper", set.REML.FALSE = TRUE, 
          reset.REML.TRUE = TRUE, log.file = file.path(tempdir(), paste("bfFixefLMER_F_log_", 
          gsub(":", "-", gsub(" ", "_", date())), ".txt", sep = ""))) 
{
  if (length(item) == 0) {
    stop("please supply a value to the ''item'' argument")
  }
  if (length(alpha) == 0) {
    stop("please supply a value to the ''alpha'' argument")
  }
  if (length(set.REML.FALSE) == 0) {
    stop("please supply a value to the ''set.REML.FALSE'' argument")
  }
  if (length(reset.REML.TRUE) == 0) {
    stop("please supply a value to the ''reset.REML.TRUE'' argument")
  }
  if (as.vector(model@call[1]) == "glmer()") {
    stop("sorry, glmer models not yet supported")
  }

  data <- model@frame
  if (llrt && set.REML.FALSE) {
    cat("setting REML to FALSE\n")
    model = update(model, . ~ ., REML = FALSE)
  }
  options(warn = 1)
  temp.dir = tempdir()
  tempdir()
  unlink(file.path(temp.dir, "temp.txt"))
  sink(file = NULL, type = "message")
  if (log.file != FALSE) 
    sink(file = log.file, split = TRUE)
  if (item != FALSE) {
    cat("checking", paste("by-", item, sep = ""), "random intercepts\n")
    model.updated <- NULL
    eval(parse(text = paste("model.updated=update(model,.~.+(1|", 
                            item, "))", sep = "")))
    if (as.vector(anova(model, model.updated)[2, "Pr(>Chisq)"]) <= 
      alpha) {
      cat("\tlog-likelihood ratio test p-value =", as.vector(anova(model, 
                                                                   model.updated)[2, "Pr(>Chisq)"]), "\n")
      cat("\tadding", paste("by-", item, sep = ""), "random intercepts to model\n")
      model = model.updated
    }
    else {
      cat("\tlog-likelihood ratio test p-value =", as.vector(anova(model, 
                                                                   model.updated)[2, "Pr(>Chisq)"]), "\n")
      cat("\tnot adding", paste("by-", item, sep = ""), 
          "random intercepts to model\n")
    }
  }
  coefs = row.names(anova(model))
  smry = pamer.fnc(model)
  smry.temp = pamer.fnc(model)
  temp = strsplit(coefs, ":")
  names(temp) = coefs
  intr.order = list()
  for (i in coefs) {
    intr.order[[i]] = length(temp[[i]])
  }
  intr.order = as.data.frame(unlist(intr.order))
  colnames(intr.order) = "Order"
  intr.order$Coef = row.names(intr.order)
  row.names(intr.order) = 1:nrow(intr.order)
  orders = sort(as.numeric(unique(intr.order$Order)), decreasing = TRUE)
  smry.temp$Order = intr.order$Order
  count = 1
  order=orders[1]
  for (order in orders) {
    cat("processing model terms of interaction level", order, 
        "\n")
    keepers = as.character(row.names(smry.temp[smry.temp$Order == 
      order, ]))
    smry.temp2 = smry.temp[keepers, ]
    if (smry.temp2[smry.temp2[, paste(p.value, ".p.val", 
        sep = "")] == max(smry.temp2[, paste(p.value, ".p.val", 
        sep = "")]), paste(p.value, ".p.val", sep = "")] <= alpha) {
      cat("\tall terms of interaction level", order, "significant\n")
    }
    else {
      while (smry.temp2[smry.temp2[, paste(p.value, ".p.val", 
        sep = "")] == max(smry.temp2[, paste(p.value, 
        ".p.val", sep = "")]), paste(p.value, ".p.val", 
        sep = "")][1] > alpha) {
        cat("\titeration", count, "\n")
        evaluated.term = row.names(smry.temp2[smry.temp2[, 
        paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, 
        paste(p.value, ".p.val", sep = "")]), ])
        cat("\t\tp-value for term", paste("\"", evaluated.term, 
         "\"", sep = ""), "=", smry.temp2[evaluated.term, 
         paste(p.value, ".p.val", sep = "")], ">",alpha, "\n")
# changed: check for occurence of the evaluated term in each higher order interaction 
#instead of checking for occurence in a combination of all terms occuring in higher order interaction                  
          hoi<- gsub("\\)","", gsub("\\(", "", intr.order[intr.order$Order > 
            order, "Coef"]))
          tt<-data.frame()
          
          for(i in 1:length(hoi)){
          tt[i,1]<-length(grep("FALSE", unique(unlist(strsplit(evaluated.term, 
             ":")) %in% unlist(strsplit(hoi[i], ":"))))) < 1}
          if (length(grep("TRUE", tt[,1])) != 0) {
#changed until here          
          cat("\t\tpart of higher-order interaction\n")
          cat("\t\tskipping term\n")
          keepers = row.names(smry.temp2)
          keepers = keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[, 
          paste(p.value, ".p.val", sep = "")] == max(smry.temp2[, 
          paste(p.value, ".p.val", sep = "")]), ])),keepers)]
          smry.temp2 = smry.temp2[keepers, ]
          smry.temp2 = na.omit(smry.temp2)
        }
        else {
          cat("\t\tnot part of higher-order interaction\n")
          m.temp <- NULL
          eval(parse(text = paste("m.temp=update(model,.~.-", 
           row.names(smry.temp2[smry.temp2[, paste(p.value, 
           ".p.val", sep = "")] == max(smry.temp2[,
           paste(p.value, ".p.val", sep = "")]),]), ")", sep = "")))
          if (llrt) {
            if (as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]) <= 
              alpha) {
              cat("\t\tlog-likelihood ratio test p-value =", 
                  as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]), 
                  "<=", alpha, "\n")
              cat("\t\tskipping term\n")
              keepers = row.names(smry.temp2)
              keepers = keepers[-grep(as.character(row.names(smry.temp2[smry.temp2[, 
              paste(p.value, ".p.val", sep = "")] == 
              max(smry.temp2[, paste(p.value, ".p.val", 
              sep = "")]), ])), keepers)]
              smry.temp2 = smry.temp2[keepers, ]
            }
            else {
              cat("\t\tlog-likelihood ratio test p.value =", 
                  as.vector(anova(model, m.temp)[2, "Pr(>Chisq)"]), 
                  ">", alpha, "\n")
             reduction<-T
            }
          }
          else {reduction=T}
          if(reduction){
            cat("\t\tremoving term\n")
            model = m.temp
            coefs = row.names(anova(model))
            if (length(coefs) == 0) {
              break
            }
            smry = pamer.fnc(model)
            smry.temp = pamer.fnc(model)
            temp = strsplit(coefs, ":")
            names(temp) = coefs
            intr.order = list()
            for (i in coefs) {
              intr.order[[i]] = length(temp[[i]])
            }
            intr.order = as.data.frame(unlist(intr.order))
            colnames(intr.order) = "Order"
            intr.order$Coef = row.names(intr.order)
            row.names(intr.order) = 1:nrow(intr.order)
            smry.temp$Order = intr.order$Order
            keepers = as.character(row.names(smry.temp[smry.temp$Order == 
              order, ]))
            smry.temp2 = smry.temp2[keepers, ]
            smry.temp2 = na.omit(smry.temp2)
          }
        }
        count = count + 1
        if (nrow(smry.temp2) == 0) {
          break
        }
      }
    }
  }
  if (llrt && reset.REML.TRUE) {
    cat("resetting REML to TRUE\n")
    model = update(model, . ~ ., REML = TRUE)
  }
  if (prune.ranefs) {
    cat("pruning random effects structure ...\n")
    split1 <- gsub(" ", "", model@call)[2]
    split2 <- unlist(strsplit(split1, "\\~"))[2]
    split2 <- gsub("\\)\\+\\(", "\\)_____\\(", split2)
    split2 <- gsub("\\(0\\+", "\\(0\\&\\&\\&", split2)
    split2 <- gsub("\\(1\\+", "\\(0\\&\\&\\&", split2)
    split3 <- unlist(strsplit(split2, "\\+"))
    split4 <- grep("\\|", split3, value = TRUE)
    split4 <- gsub("\\&\\&\\&", "\\+", split4)
    split5 <- unlist(strsplit(split4, "_____"))
    m.ranefs <- split5
    coefs = row.names(anova(model))
    m.ranef.variables <- gsub("\\((.*)\\|.*", "\\1", m.ranefs)
    m.ranef.variables <- gsub(".\\+(.*)", "\\1", m.ranef.variables)
    m.ranef.variables <- m.ranef.variables[m.ranef.variables != 
      "1"]
    m.ranef.variables <- m.ranef.variables[m.ranef.variables != 
      "0"]
    ranef.to.remove <- vector("character")
    for (m.ranef.variable in m.ranef.variables) {
      if (!m.ranef.variable %in% coefs) {
        ranef.to.remove <- c(ranef.to.remove, m.ranef.variable)
      }
    }
    if (length(ranef.to.remove) > 0) {
      rtr <- vector("character")
      for (iii in ranef.to.remove) {
        rtr <- c(rtr, grep(iii, m.ranefs, value = TRUE))
        cat(paste("\t", iii, " in random effects structure, \n\tbut not in fixed effects structure\n", 
                  sep = ""))
        cat("\t\tremoving", iii, "from model ...\n")
      }
      eval(parse(text = paste("model<-update(model,.~.-", 
                              paste(rtr, collapse = "-"), ",data=data)", sep = "")))
    }
    else {
      cat("\tnothing to prune\n")
    }
  }
  options(warn = 0)
  sink(file = NULL, type = "message")
  if (log.file != FALSE) {
    cat("log file is", log.file, "\n")
    sink(file = NULL)
  }
  return(model = model)
}
