tidy.gls <- function(x, conf.int = FALSE, conf.level = .95,
                    exponentiate = FALSE, quick = FALSE, ...) {
    if (quick) {
        co <- stats::coef(x)
        ret <- data.frame(term = names(co), estimate = unname(co))
        return(process_lm(ret, x, conf.int = FALSE, exponentiate = exponentiate))
    }
    s <- summary(x)
    ret <- tidy.summary.lm(s)
    
    process_lm(ret, x, conf.int = conf.int, conf.level = conf.level,
         exponentiate = exponentiate)
}


#' @rdname lm_tidiers
#' @export
tidy.summary.lm <- function(x, ...) {
    co <- stats::coef(x)
    nn <- c("estimate", "std.error", "statistic", "p.value")
    if (inherits(co, "listof")) {
        # multiple response variables
        ret <- plyr::ldply(co, fix_data_frame, nn[1:ncol(co[[1]])],
                           .id = "response")
        ret$response <- stringr::str_replace(ret$response, "Response ", "")
    } else {
        ret <- fix_data_frame(co, nn[1:ncol(co)])
    }
    
    ret
}

process_lm <- function(ret, x, conf.int = FALSE, conf.level = .95,
                       exponentiate = FALSE) {
    if (exponentiate) {
        # save transformation function for use on confidence interval
        if (is.null(x$family) ||
            (x$family$link != "logit" && x$family$link != "log")) {
            warning(paste("Exponentiating coefficients, but model did not use",
                          "a log or logit link function"))
        }
        trans <- exp
    } else {
        trans <- identity
    }
    
    if (conf.int) {
        # avoid "Waiting for profiling to be done..." message
        CI <- suppressMessages(stats::confint(x, level = conf.level))
        colnames(CI) = c("conf.low", "conf.high")
        ret <- cbind(ret, trans(unrowname(CI)))
    }
    ret$estimate <- trans(ret$estimate)
    
    ret
}
