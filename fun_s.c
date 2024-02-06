/*****************
 AC - OpenMP -- SERIE
 fun_s.c
 rutinas que se utilizan en el modulo grupopal_s.c
******************/
#include <math.h>
#include <stdlib.h>
#include "defineg.h" // definiciones

/***********************
 1 - Funcion para calcular la distancia euclidea entre dos vectores
 Entrada: 2 elementos con NDIM caracteristicas (por referencia)
 Salida:  distancia (double)
************************/
double gendist (float *vec1, float *vec2){
	double sumaCuadrados = 0.00;
	for(int i = 0; i < NDIM; i++)
		sumaCuadrados += pow((vec1[i] - vec2[i]), 2);
	return sqrt(sumaCuadrados);
}

/*****************************
 2 - Funcion para calcular el grupo (cluster) mas cercano (centroide mas cercano)
 Entrada: nvec  numero de vectores, int
          mvec  vectores, una matriz de tamanno MAXV x NDIM, por referencia
          cent  centroides, una matriz de tamanno ngrupos x NDIM, por referencia
 Salida:  popul grupo mas cercano a cada elemento, vector de tamanno MAXV, por ref.
****************************/
void grupo_cercano (int nvec, float mvec[][NDIM], float cent[][NDIM],
		int *popul){
	int centroide;
	float distancia_minima, distancia_actual;
	for(int i = 0; i < nvec; i++){
		distancia_minima = 10.00;
		centroide = -1;
		for(int j = 0; j < ngrupos; j++){
			distancia_actual = gendist(mvec[i], cent[j]);
			if(distancia_actual <= distancia_minima){
				centroide = j;
				distancia_minima = distancia_actual; 
			}
		}
		popul[i] = centroide;
	}
}

/*****************************
 3 - Funcion para calcular la calidad de la particion de clusteres.
     Ratio entre a y b.
     El termino a corresponde a la distancia intra-cluster.
     El termino b corresponde a la distancia inter-cluster.
 Entrada: mvec    vectores, una matriz de tamanno MAXV x NDIM, por referencia
          listag  vector de ngrupos structs (informacion de grupos generados), por ref.
          cent    centroides, una matriz de tamanno ngrupos x NDIM, por referencia
 Salida:  valor del CVI (double): calidad/bondad de la particion de clusters
******************************/
double silhouette_simple(float mvec[][NDIM], struct lista_grupos *listag, float cent[][NDIM], float a[]) {
    double distInter, distIntra, sumaDistancias, sumaTotal, sumaRatios = 0;
	int i, j, k, p, *vecgAux, nvecgAux;

	for (k = 0; k < ngrupos; k++) {

		nvecgAux = listag[k].nvecg;
		vecgAux = listag[k].vecg;

		sumaDistancias = 0;
		distIntra = 0;

		if (nvecgAux > 1) {
				for(i = 0; i < nvecgAux; i++)
					for(j = 0; j < nvecgAux; j++ )
						sumaDistancias += gendist (mvec[vecgAux[i]], mvec[vecgAux[j]]);

			distIntra = sumaDistancias / (nvecgAux * (nvecgAux-1));
		} a[k] = distIntra;

		sumaDistancias = 0;
			for(p = 0; p < ngrupos; p++)
				sumaDistancias += gendist(cent[k], cent[p]);
			
		distInter = sumaDistancias / (ngrupos - 1);
		
		sumaRatios += (distInter - distIntra) / fmax(distInter, distIntra);
	} return sumaRatios / ngrupos;
}
/********************************
 4 - Funcion para realizar el analisis de campos UNESCO
 Entrada:  listag   vector de ngrupos structs (informacion de grupos generados), por ref.
           mcam     campos, una matriz de tamaño MAXV x NCAM, por referencia
 Salida:   info_cam vector de NCAM structs (informacion del analisis realizado), por ref.
*******************************/
float medianaOrdenada(float *vec, int num) {
    double aux;
    for (int i = 0; i < num; i++) {
        for (int j = 0; j < (num-i-1); j++) {
            if (vec[j] < vec[j+1]) {
                aux = vec[j];
                vec[j] = vec[j+1];
                vec[j+1] = aux;
            }
        }
    }
    if (num % 2 == 0) {
		return fmax(vec[num/2], vec[(num/2)-1]);
    } else  return vec[num/2];
}

void analisis_campos (struct lista_grupos *listag, float mcam[][NCAM], struct analisis *info_cam) {
	// Realizar el analisis de campos UNESCO en los grupos:
	//    mediana maxima y el grupo en el que se da este maximo (para cada campo)
	//    mediana minima y su grupo en el que se da este minimo (para cada campo)
	int nvecg, gcmin, gcmax;
	float *vecGrupo;
	double min, max, mediana;
	
	for (int c = 0; c < NCAM; c++) {
		min = 1;
		max = 0;
		for (int k = 0; k < ngrupos; k++) {
			nvecg = listag[k].nvecg;
			if(nvecg > 0) {
				vecGrupo = calloc(nvecg, sizeof(float));
				for (int i = 0; i < nvecg; i++)
					vecGrupo[i] = mcam[listag[k].vecg[i]][c];
				mediana = medianaOrdenada(vecGrupo, nvecg);
				if (mediana < min) {
					gcmin = k;
					min = mediana;
				}
				if (mediana > max) {
					gcmax = k;
					max = mediana;
				}
				free(vecGrupo);
			}
		}
		info_cam[c].mmin = min;
		info_cam[c].mmax = max;
		info_cam[c].gmin = gcmin;
		info_cam[c].gmax = gcmax;
	}
}
/*************
   OTRAS FUNCIONES DE LA APLICACION
**************/
void inicializar_centroides(float cent[][NDIM]){
	int i, j;
	float rand_val;
	srand (147);
	for (i=0; i<ngrupos; i++)
		for (j=0; j<NDIM/2; j++){
			rand_val = ((rand() % 10000) / 10000.0)*2-1;
			cent[i][j] = rand_val;
			cent[i][j+(NDIM/2)] = cent[i][j];
		}
}

int nuevos_centroides(float mvec[][NDIM], float cent[][NDIM], int popul[], int nvec){
	int i, j, fin;
	double discent;
	double additions[ngrupos][NDIM+1];
	float newcent[ngrupos][NDIM];

	for (i=0; i<ngrupos; i++)
		for (j=0; j<NDIM+1; j++)
			additions[i][j] = 0.0;

	// acumular los valores de cada caracteristica; numero de elementos al final
	for (i=0; i<nvec; i++){
		for (j=0; j<NDIM; j++) additions[popul[i]][j] += mvec[i][j];
		additions[popul[i]][NDIM]++;
	}

	// calcular los nuevos centroides y decidir si el proceso ha finalizado o no (en funcion de DELTA)
	fin = 1;
	for (i=0; i<ngrupos; i++){
		if (additions[i][NDIM] > 0) { // ese grupo (cluster) no esta vacio
			// media de cada caracteristica
			for (j=0; j<NDIM; j++)
				newcent[i][j] = (float)(additions[i][j] / additions[i][NDIM]);

			// decidir si el proceso ha finalizado
			discent = gendist (&newcent[i][0], &cent[i][0]);
			if (discent > DELTA1){
				fin = 0;  // en alguna centroide hay cambios; continuar
			}

			// copiar los nuevos centroides
			for (j=0; j<NDIM; j++)
				cent[i][j] = newcent[i][j];
		}
	}
	return fin;
}