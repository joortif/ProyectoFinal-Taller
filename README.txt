#PROYECTO FINAL TALLER TRANSVERSAL I: PROGRAMACIÓN Y PROCESO DE INFORMACIÓN
#TRADUCCIÓN DEL PROYECTO INTIAL ORBIT DETERMINATION A C++

Autor: Joaquín Ortiz de Murua Ferrero

El presente archivo es el README para el proyecto de traducción de MATLAB a C++
del programa para la determinación inicial de órbitas "Initial Orbit Determination".

ESTRUCTURA DEL PROGRAMA

El programa está estructurado en 4 carpetas principales:
-include: Contiene las cabeceras de las funciones utilizadas en el proyecto (*.h).
-data: Contiene los archivos de texto plano (*.txt) necesarios para obtener los valores de algunas variables y componentes utilizados.
-src: Contiene los archivos de código fuente de C++ utilizados donde se definen las funciones del proyecto principal (archivos *.cpp).
-html: Contiene la documentación del proyecto generada automáticamente con la herramienta Doxygen y los gráficos sobre la ejecución de las funciones generado por Graphviz.

Los ficheros adicionales que se encuentran en la carpeta son:
-EKF_GEOS3.cpp: Archivo con el programa principal del proyecto
-EKF_Tests.cpp: Archivo de pruebas unitarias realizadas sobre las funciones utilizadas en el proyecto.

INSTRUCCIONES DE COMPILACIÓN Y EJECUCIÓN

1. Descargar la carpeta con los archivos del proyecto y descomprimirla en el sistema local.
2. Asegurarse de que se dispone de un compilador C++ GCC compatible con la versión 9.2.0 y de las bibliotecas cmath, valarray, fstream e iostream de C++ instaladas en el sistema.
3. Abrir una terminal del sistema operativo que se esté utilizando y navegar hasta el directorio raíz del proyecto.
4. Ejecutar el comando de compilación siguiente para el programa principal:
	
		g++ -o orbitDetermination EKF_GEOS3.cpp src/*.cpp 

   Si se desea ejecutar el archivo de pruebas el comando a ejecutar es el siguiente:

    	g++ -o orbitDeterminationTests EKF_Tests.cpp src/*.cpp 

5. Una vez compilado con éxito, ejecutar el programa con el siguiente comando:

	./orbitDetermination
