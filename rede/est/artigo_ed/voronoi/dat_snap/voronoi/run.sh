#!/bin/bash
#
#run.sh - Roda o programa vor.c para diferentes phi-X
#
# Criado em     : 2024/03/07
# Atualizado		: 2024/03/07
# Autor		      : Matheus <matheuslucasgp@gmail.com>
#----------------------------------------
# Notas:
#
#----------------------------------------
# Exemplo:
# ./run.sh 15 50
# 	Roda o programa vor.out com seed 15 até seed 50
#----------------------------------------
# Histórico:
# 	V1.0 - Roda o programa fazendo plot
#----------------------------------------

#---Variáveis---
SDI=$1
SDF=$2

#---Verificação--
if [ $# -ne 2 ]; then
	echo "Uso: $0 <Sd I> <Sd F>"
	echo "<Sd I>: Seed inicial"
	echo "<Sd F>: Seed final"
  exit 1
fi

for (( I = SDI; I <= SDF ; I++ )); do
	echo -ne "Executando plot $I/$SDF\r"	
	make
	./vor.out $I > plt/phi-$I.png
done
