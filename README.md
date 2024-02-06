Este proyecto se centra en la paralelización de una aplicación que realiza clustering de palabras basadas en embeddings y análisis de su relación con campos definidos por la UNESCO. Inspirada en el procesamiento del lenguaje natural y el aprendizaje automático, la aplicación consta de dos fases principales.

En la Fase 1, la aplicación lee dos archivos: uno que contiene vectores de palabras y otro que asocia estos vectores con campos UNESCO. Luego, genera centroides iniciales aleatorios para K grupos. Utilizando cálculos de distancias euclidianas, asigna cada palabra al grupo más cercano y actualiza los centroides. Este proceso se repite hasta converger o alcanzar un límite de iteraciones.

En la Fase 2, se analiza la relación de los grupos de palabras con los campos UNESCO. Se determina qué grupo está más cerca y cuál está más lejos de cada campo, utilizando la mediana de las distancias. Estas distancias se ordenan de menor a mayor utilizando un algoritmo simple como el de la burbuja, y luego se calcula la mediana para cada campo.

La paralelización del proyecto se ha llevado a cabo mediante OpenMP para aprovechar la capacidad de cómputo de sistemas multiprocesador. Se han identificado y paralelizado las secciones críticas de la aplicación, mejorando así la eficiencia del proceso de clustering y análisis de campos.
