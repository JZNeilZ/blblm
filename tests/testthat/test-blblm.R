test_that("Output of blblm is S3 class object and
          the shape of coefficient is correct", {
            fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 2, B = 100)
            co <- coef(fit)
            ci <- confint(fit)

            expect_s3_class(fit, "blblm")
            expect_equal(length(co), 4)
            expect_equal(dim(ci), c(3,2))

            fit1 <- blbglm(mpg ~ wt + hp, data = mtcars, family = "gaussian", m = 2, B = 100)
            co1 <- coef(fit1)
            ci1 <- confint(fit1)

            expect_s3_class(fit1, "blblm")
            expect_equal(length(co1), 3)
            expect_equal(dim(ci1), c(2,2))
})
