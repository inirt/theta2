
install.packages('./', type = 'source', repos = NULL)

library(instrument)

data(familyrisk)

model = 'theta1 = c(3:16)
         theta2 = c(3:16)
         theta3 = c(3:16)'

fit =
  instrument::instrument(
    data = familyrisk,
    model = model,
    iter_warmup = 10,
    iter_sampling = 20,
    seed = 12345
    )





# convenient install using url!
# install.packages("https://github.com/hadley/devtools/archive/v1.7.0.tar.gz",
#                  repos=NULL, method="libcurl")
