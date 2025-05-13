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

setVariable(matlab, X = Meat$x)
#setVariable(matlab, X = dat)

evaluate(matlab, "whos")

UFS_Methods <- list(
  "CFS" = FALSE,
  "Laplacian" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  "RUFS" = TRUE
)

UFS_Methods_Slow <- list(
  "MCFS" = TRUE,
  "LLCFS" = TRUE,
  "FSASL" = TRUE,
  "SOCFS" = TRUE
)

for (UFS_Method in names(UFS_Methods))  {
  print(UFS_Method)
  send_input <- UFS_Methods[[UFS_Method]]
  cmd <- sprintf('Result = Auto_UFSTool(X, "%s");', UFS_Method)
  evaluate(matlab, "rng(42);")
  # Execute AutoHotkey to send option 1 keystroke (default values) to MATLAB console
  if (send_input) {
    system('cmd /c start /B "C:/Program Files/AutoHotkey/v2/AutoHotkey.exe" "C:/Users/anton/OneDrive/Documentos/send_input.ahk"')
  }
  # Ahora ejecutar la función en MATLAB
  evaluate(matlab, cmd)
  result <- getVariable(matlab, "Result")
  print(str(result))
}

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


# Iniciar servidores en distintos puertos
Matlab$startServer(port = 9999)
Matlab$startServer(port = 10000)

# Crear objetos MATLAB conectados a cada servidor
matlab1 <- Matlab(port = 9999)
matlab2 <- Matlab(port = 10000)

# Conectar ambos objetos a los servidores
open(matlab1)
open(matlab2)

# Ver información de cada conexión
print(matlab1)
print(matlab2)

evaluate(matlab1, "set(groot, 'Name', 'MATLAB - Servidor 9999');")

close(matlab1)
close(matlab2)



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

evaluate(matlab, "disp('1');") # Enviar '1' a la consola de MATLAB


