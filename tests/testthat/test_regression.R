library(grpack)
data(gen_01)
lm1 <- reg(y~x_1-1, data = gen_01)
lm2 <- reg(y~x_1, data = gen_01)
lm3 <- reg(y~x_1, data = gen_01, weights = weights)
lm4 <- reg(y~x_1+x_2, data = gen_01)
lm5 <- reg(y~x_1+x_2, data = gen_01, weights = weights)


context("Testing flat regression")


test_that("Test regression est on generate data gen_dat_01.rda", {
  expect_is(lm1, "reg")
  expect_is(lm2, "reg")
  expect_is(lm3, "reg")
  expect_is(lm4, "reg")
  expect_is(lm5, "reg")

  expect_equal(coef(lm1), c(x_1=-0.00545594686866))
  expect_equal(coef(lm2), structure(c(0.108852115808104, -0.000932369740880879),
                                    .Names = c("(Intercept)", "x_1")))
  expect_equal(coef(lm3), structure(c(0.21247598889688, -0.0901806422487909),
                                    .Names = c("(Intercept)", "x_1")))

  expect_equal(coef(lm5), structure(c(0.168349625509692, -0.0933352340816591,
                                      0.0873893119812259),
                                    .Names = c("(Intercept)", "x_1", "x_2")))
})

test_that("Test regression vcov on generate data gen_dat_01.rda", {
  expect_equal(sqrt(diag(vcov(lm5, type = "const"))),
               structure(c(0.172868824202954, 0.0935384796984826,
                           0.292468328116345),
                         .Names = c("(Intercept)", "x_1", "x_2")))

  expect_equal(sqrt(diag(vcov(lm5, type = "HC1"))),
    structure(c(0.18234976476785, 0.0860612203065056,
                0.340011160549175),
            .Names = c("(Intercept)", "x_1", "x_2")))

  expect_equal(sqrt(diag(vcov(lm5, type = "HC2"))),
               structure(c(0.183626662423246, 0.0870453887156868,
                           0.342282200513398),
                         .Names = c("(Intercept)", "x_1", "x_2")))

  expect_equal(sqrt(diag(vcov(lm5, type = "HC3"))),
    structure(c(0.187777796737193, 0.0894191619373881,
                0.349906505986469),
               .Names = c("(Intercept)", "x_1", "x_2")))

})

