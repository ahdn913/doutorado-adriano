#!/bin/bash
#
#vis.sh - Shell para visualizar diferentes redes pelo gnuplot
#
# Criado em     : 2023/11/08
# Atualizado		: 2023/11/08
# Autor		      : Matheus <matheuslucasgp@gmail.com>
#----------------------------------------
# Notas:
#
#----------------------------------------
# Exemplo:
# ./vis.sh 1 8  gera todos os ind. ativos de rede 8x8 triangular 
#----------------------------------------
# Histórico:
#
#----------------------------------------

#---CHAVES---
if [[ $1 -eq 1 ]]; then TRIAN=1; fi
if [[ $1 -eq 2 ]]; then HONEY=1; fi
if [[ $1 -eq 3 ]]; then KAGOM=1; fi
if [[ $1 -eq 4 ]]; then 					
	TRIAN=1;
 	HONEY=1;
 	KAGOM=1;
fi

#---VERIFICAÇÕES---
#quantidade de termos
if [[ $# -ne 2 ]]; then 
	echo "Uso: $0 <Tipo rede> <Tamanho rede>"
	echo "Tipo = 1 : Rede Triangular"
	echo "Tipo = 2 : Rede Favo de mel"
	echo "Tipo = 3 : Rede Kagomé"
	echo "Tipo = 4 : Todas as redes" 
	exit 1
fi 

#---Rede Triangular---

if [[ $TRIAN -eq 1 ]]; then
	echo "Executando rede triangular $2 X $2..."
	make triang
	for (( i = 0; i < $(($2*$2)); i++ ))
	do				
		./rede-triang.out $2 $i
		cd plt/
		gnuplot -c rede-triang.plt $2 $i
		cd ..
	done
fi

#---Rede Honey---

if [[ $HONEY -eq 1 ]];
then
	echo "Executando rede honeycomb $2 X $2..."
	make honey
	for (( i = 0; i < $(($2*$2)); i++ ))
	do				
		./rede-honeycomb.out $2 $i
		cd plt/
		gnuplot -c rede-honey.plt $2 $i
		cd ..
	done
fi

#---Rede Kagomé---

if [[ $KAGOM -eq 1 ]];
then
	echo "Executando rede kagomé $2 X $2..."
	make kagome
	for (( i = 0; i < $((3*$2*$2/4)); i++ ))
	do				
		./rede-kagome.out $2 $i
		cd plt/
		gnuplot -c rede-kagome.plt $2 $i
		cd ..
	done
fi
