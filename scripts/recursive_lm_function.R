recursive_lm = function(outcome, predictor, data, formula){
  # 1) outcome list of strings with the outcomes (Ys)
  # 2) predictor string with the predictor whose beta coefficient would be obtained (Xs)
  # 3) data dataframe with outcomes, predictors and element of formula (data)
  # 4) formula string with the formula string for lm function 
  #### lm object with standardized beta coefficient 
  object1 = lm.beta(lm(as.formula(paste(outcome,"~", 
                                        predictor, formula)), 
                       data = data))
  object2 = lm(as.formula(paste(outcome,"~", 
                                predictor, formula)), 
               data = data)
  #### list in which each element is a summary dataframe about the lm model
  tidy_object1 =  broom::tidy(object1) %>%
    dplyr::select(term, std_estimate) %>%
    dplyr::rename(predictor = term) %>%
    dplyr::rename(std_coefficient = std_estimate)
  tidy_object2 = broom::tidy(object2, conf.int = TRUE) %>% 
    dplyr::filter(term == predictor) %>% 
    dplyr::rename(predictor = term) %>% 
    dplyr::mutate(outcome = outcome, 
                  formula = paste(attr(object2$terms,"term.labels"), collapse = " ")) %>% 
    dplyr::relocate(outcome, .before = estimate) %>% 
    dplyr::rename(beta_coefficient = estimate, 
                  # std_coefficient = std_estimate, 
                  std_error = std.error,
                  t_statistic = statistic,
                  CI_low = conf.low, 
                  CI_high = conf.high)
  tidy_object = merge(tidy_object1, tidy_object2, by = "predictor") %>% 
    dplyr::relocate(outcome, .before = predictor) %>% 
    dplyr::relocate(beta_coefficient, .before = std_coefficient)
  return(tidy_object)
}