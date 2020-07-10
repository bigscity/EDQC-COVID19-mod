my_data <- read.csv("./simuresult/resultsimu111.csv")
res <- cor.test(my_data$reported_I, my_data$simulation_I, method='pearson')
print(res)