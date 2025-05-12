# Connection to MATLAB

library(R.matlab)
library(diceRplus)
library(RPEClust)

Matlab$startServer()
matlab <- Matlab()
print(matlab)

isOpen <- open(matlab)
if (!isOpen)
  print("MATLAB server is not running: waited 30 seconds.")
print(matlab)

data(meats)

str(meats$class)

dat <-  meats[, !names(meats) %in% "class"]

ref.cl <- as.integer(meats$class)

# setVariable(matlab, X = Meat$x)
setVariable(matlab, X = dat)

evaluate(matlab, "whos")

evaluate(matlab, "rng(42);")
evaluate(matlab, "Result = Auto_UFSTool(X, 'RUFS');")
result <- getVariable(matlab, "Result")
# Extraer el vector numérico
ranking <- result$Result[[1]]
# Aplana el array:
ranking <- as.vector(ranking)

# Ahora puedes ver los primeros elementos con:
head(ranking)

# O imprimir todo el ranking:
print(ranking)

close(matlab)


######




# Carga el archivo COIL20.mat
evaluate(matlab, "load('COIL20.mat');")
# Ejecuta Auto_UFSTool
evaluate(matlab, "rng(42);")
evaluate(matlab,"Result = Auto_UFSTool(X, 'RUFS');")


result <- getVariable(matlab, "Result")
# Extraer el vector numérico
ranking <- result$Result[[1]]
# Aplana el array:
ranking <- as.vector(ranking)

# Ahora puedes ver los primeros elementos con:
head(ranking)

# O imprimir todo el ranking:
print(ranking)

# String to send option 1 to MATLAB Server
codigo <- paste(
  "robot = java.awt.Robot;",
  "pause(4);",
  "robot.keyPress(java.awt.event.KeyEvent.VK_1);",
  "robot.keyRelease(java.awt.event.KeyEvent.VK_1);",
  "robot.keyPress(java.awt.event.KeyEvent.VK_ENTER);",
  "robot.keyRelease(java.awt.event.KeyEvent.VK_ENTER);",
  "Result = Auto_UFSTool(X, 'RUFS');",
  sep = " "
)

path_MATLAB = "C:/Users/anton/OneDrive/Documentos/MATLAB"

fileConn <- file(paste0(path_MATLAB, "/AutoRun.m"))
writeLines(codigo, fileConn)
close(fileConn)





evaluate(matlab, "
options = struct('MaxIter', 10, 'epsilon', 1e-4, 'nu', 1, 'alpha', 1, 'beta', 1);
Result = Auto_UFSTool(X, 'RUFS', options);
")

evaluate(matlab, "Result;")
result <- getVariable(matlab, "Result")
close(matlab)


evaluate(matlab, "A = 1+2;", "B = ones(2, 20);")
evaluate(matlab, "A")
data <- getVariable(matlab, c("A", "B"))
cat("Received variables:\n")
str(data)
close(matlab)
