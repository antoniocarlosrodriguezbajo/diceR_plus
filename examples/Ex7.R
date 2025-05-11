# Connection to MATLAB

library(R.matlab)
Matlab$startServer()

matlab <- Matlab()

print(matlab)

isOpen <- open(matlab)
if (!isOpen)
  throw("MATLAB server is not running: waited 30 seconds.")
print(matlab)

# Carga el archivo COIL20.mat
evaluate(matlab, "load('COIL20.mat');")
# Ejecuta Auto_UFSTool
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
