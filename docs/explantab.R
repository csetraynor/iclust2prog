explanationrow <- c("Accuracy of model prediction", "perfect", "strong", "substantial", "moderate", "borderline")
explanation <- data.frame(BS =  "Accuracy of model prediction",
                          '0' =  "perfect",
                          "0-0.025" = "overwhelming",
                          "0.025 - 0.05" = "strong",
                          "0.05 - 0.1" = "substantial",
                          "0.1 - 0.2" = "moderate",
                          "0.2 - 0.25" = "borderline")
require(stargazer)
stargazer(explanationdata  , type = "latex", summary = FALSE, digits.extra = 3,
          digits = 3, digit.separator = ".",
          rownames = F)


explanationdata<- data.frame("Accuracy of model prediction" = c( "perfect", "overwhelming", "strong","substantial", "moderate", "borderline", "trivial"),
                             BS = c(0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25))
