#include "prepbox.cpp"
#include <stdlib.h>
#include <stdio.h>

#define RAIO_ATUACAO 1.5 // Em relacao ao intruso
#define PERIODO 1000

struct CELULA{
	unsigned int contador;
    int *lista;
};

class PREDITOR_CORRETOR : public CAIXA {
    public:
        double dt, gravidade, p0y[4], gn;
    public:
        unsigned long long int passos;
        PREDITOR_CORRETOR();
        PREDITOR_CORRETOR(double, double, unsigned int);
        PREDITOR_CORRETOR(int, int, double**);
        PREDITOR_CORRETOR(int, int, double**, int, int**, double*);
        ~PREDITOR_CORRETOR();
        int Velocidade_Aleatoria();
        int Procura_Vizinhos(unsigned int);
		int Procura_Vizinhos2(unsigned int);
        int Detectar_Contatos(unsigned int);
        int Verifica_Equilibrio();
        double Coordenacao(int);
        double Interpenetracao();
        double Compactacao();
        double Pressao_Parede(int);
        double Energia_Cinetica();
        double Energia_Rotacional();
        double Energia_Elastica();
        double Energia_Gravitacional();
        double Energia_Cinetica(int);
        double Energia_Rotacional(int);
        double Energia_Elastica(int);
        double Energia_Gravitacional(int);
        int Preditor(unsigned int);
        int Subir_Parede(double);
        int Vibra_Parede(double, double);
        int Calculo_Forca(unsigned int);
        int Corretor(unsigned int);
        int Salvar_Configuracao();
        int Salvar_Configuracao(const char *);
        double Area_Interpenetracao(int);
        double Area_Interpenetracao_Total();
};

PREDITOR_CORRETOR::PREDITOR_CORRETOR(int p, int n, double **dados, int numero_contatos, int **contatos, double *f_contatos) : CAIXA(p, n, dados){
    passos = 0;
    double  menor_massa = Intruso.massa, maior_mola = Intruso.mola.normal, menor_mpk = Intruso.massa/Intruso.mola.normal, menor_mvk = Intruso.massa*Intruso.mola.normal;

    for(int i = 0; i < 4; i++){
        menor_mpk = (menor_mpk > Bordas[i].massa/Bordas[i].mola.normal) ? Bordas[i].massa/Bordas[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Bordas[i].massa*Bordas[i].mola.normal) ? Bordas[i].massa*Bordas[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Bordas[i].massa) ? Bordas[i].massa : menor_massa;
        maior_mola = (maior_mola < Bordas[i].mola.normal) ? Bordas[i].mola.normal : maior_mola;
        p0y[i] = Bordas[i].Linear.Posicao.y;
    }

    for(int i = 0; i < numero_graos; i++){
        menor_mpk = (menor_mpk > Graos[i].massa/Graos[i].mola.normal) ? Graos[i].massa/Graos[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Graos[i].massa*Graos[i].mola.normal) ? Graos[i].massa*Graos[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Graos[i].massa) ? Graos[i].massa : menor_massa;
        maior_mola = (maior_mola < Graos[i].mola.normal) ? Graos[i].mola.normal : maior_mola;
    }

    dt = sqrt(menor_mpk) /10;
    gn = 2*sqrt(menor_mvk);
    gravidade = 1;
}

PREDITOR_CORRETOR::PREDITOR_CORRETOR(int p, int n, double **dados) : CAIXA(p, n, dados){
    passos = 0;
    double  menor_massa = Intruso.massa, maior_mola = Intruso.mola.normal, menor_mpk = Intruso.massa/Intruso.mola.normal, menor_mvk = Intruso.massa*Intruso.mola.normal;

    for(int i = 0; i < 4; i++){
        menor_mpk = (menor_mpk > Bordas[i].massa/Bordas[i].mola.normal) ? Bordas[i].massa/Bordas[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Bordas[i].massa*Bordas[i].mola.normal) ? Bordas[i].massa*Bordas[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Bordas[i].massa) ? Bordas[i].massa : menor_massa;
        maior_mola = (maior_mola < Bordas[i].mola.normal) ? Bordas[i].mola.normal : maior_mola;
        p0y[i] = Bordas[i].Linear.Posicao.y;
    }

    for(int i = 0; i < numero_graos; i++){
        menor_mpk = (menor_mpk > Graos[i].massa/Graos[i].mola.normal) ? Graos[i].massa/Graos[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Graos[i].massa*Graos[i].mola.normal) ? Graos[i].massa*Graos[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Graos[i].massa) ? Graos[i].massa : menor_massa;
        maior_mola = (maior_mola < Graos[i].mola.normal) ? Graos[i].mola.normal : maior_mola;
    }

    dt = sqrt(menor_mpk) /10;
    gn = 2*sqrt(menor_mvk);
    gravidade = 1;
}

PREDITOR_CORRETOR::PREDITOR_CORRETOR() : CAIXA(){
    passos = 0;
    double  menor_massa = Intruso.massa, maior_mola = Intruso.mola.normal, menor_mpk = Intruso.massa/Intruso.mola.normal, menor_mvk = Intruso.massa*Intruso.mola.normal;

    for(int i = 0; i < 4; i++){
        menor_mpk = (menor_mpk > Bordas[i].massa/Bordas[i].mola.normal) ? Bordas[i].massa/Bordas[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Bordas[i].massa*Bordas[i].mola.normal) ? Bordas[i].massa*Bordas[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Bordas[i].massa) ? Bordas[i].massa : menor_massa;
        maior_mola = (maior_mola < Bordas[i].mola.normal) ? Bordas[i].mola.normal : maior_mola;
        p0y[i] = Bordas[i].Linear.Posicao.y;
    }

    for(int i = 0; i < numero_graos; i++){
        menor_mpk = (menor_mpk > Graos[i].massa/Graos[i].mola.normal) ? Graos[i].massa/Graos[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Graos[i].massa*Graos[i].mola.normal) ? Graos[i].massa*Graos[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Graos[i].massa) ? Graos[i].massa : menor_massa;
        maior_mola = (maior_mola < Graos[i].mola.normal) ? Graos[i].mola.normal : maior_mola;
    }

    dt = sqrt(menor_mpk) /10;
    gn = 2*sqrt(menor_mvk);
    gravidade = 1;
}

PREDITOR_CORRETOR::PREDITOR_CORRETOR(double l, double p, unsigned int n) : CAIXA(l, p, n){
    passos = 0;
    double  menor_massa = Intruso.massa, maior_mola = Intruso.mola.normal, menor_mpk = Intruso.massa/Intruso.mola.normal, menor_mvk = Intruso.massa*Intruso.mola.normal;

    for(int i = 0; i < 4; i++){
        menor_mpk = (menor_mpk > Bordas[i].massa/Bordas[i].mola.normal) ? Bordas[i].massa/Bordas[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Bordas[i].massa*Bordas[i].mola.normal) ? Bordas[i].massa*Bordas[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Bordas[i].massa) ? Bordas[i].massa : menor_massa;
        maior_mola = (maior_mola < Bordas[i].mola.normal) ? Bordas[i].mola.normal : maior_mola;
        p0y[i] = Bordas[i].Linear.Posicao.y;
    }

    for(int i = 0; i < numero_graos; i++){
        menor_mpk = (menor_mpk > Graos[i].massa/Graos[i].mola.normal) ? Graos[i].massa/Graos[i].mola.normal : menor_mpk;
        menor_mvk = (menor_mvk > Graos[i].massa*Graos[i].mola.normal) ? Graos[i].massa*Graos[i].mola.normal : menor_mvk;
        menor_massa = (menor_massa > Graos[i].massa) ? Graos[i].massa : menor_massa;
        maior_mola = (maior_mola < Graos[i].mola.normal) ? Graos[i].mola.normal : maior_mola;
    }

    dt = sqrt(menor_mpk) /10;
    gn = 2*sqrt(menor_mvk);
    gravidade = 1;
    Velocidade_Aleatoria();
}

PREDITOR_CORRETOR::~PREDITOR_CORRETOR(){
}

int PREDITOR_CORRETOR::Velocidade_Aleatoria(){
     double angulo = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for ordered
#endif
    for(int i = 0; i < numero_graos; i++){
#ifdef _OPENMP
        #pragma omp ordered
#endif
        {
            angulo = 2.0 * PI * random(&semente);
            Graos[i].Linear.Velocidade = {random(&semente) * Graos[i].Geometria.raio * cos(angulo), random(&semente) * Graos[i].Geometria.raio * sin(angulo)};
        }
    }

    return 0;
}

int PREDITOR_CORRETOR::Procura_Vizinhos2(unsigned int CONTROLE){
    Intruso.vizinhos[0] = 0;

    double minx = Bordas[0].Linear.Posicao.x +Bordas[0].Geometria.raio;
    double miny = Bordas[0].Linear.Posicao.y +Bordas[0].Geometria.raio;
    double maxx = Bordas[3].Linear.Posicao.x -Bordas[3].Geometria.raio;
    double maxy = Bordas[3].Linear.Posicao.y -Bordas[3].Geometria.raio;

    double raio_maior = 0.0;
    for (int i = 0; i < numero_graos; i ++){
        raio_maior = (Graos[i].Geometria.raio > raio_maior) ? Graos[i].Geometria.raio : raio_maior;
    }
//    raio_maior = (Intruso.Geometria.raio > 3.0*raio_maior) ? Intruso.Geometria.raio/1.5 : raio_maior;

    double L = maxx-minx;
    double M = maxy-miny;
    double l = raio_maior*3; // Largura da caixa de busca do procura de vizinhos 2
    double m = l;             // Comprimento da caixa de busca do procura de vizinhos 2
    int Ll = int(L/l) +1, Mm = int(M/m) +1;
    CELULA Celula[Ll][Mm];

    for(int i = 0; i < Ll; i++)
        for(int j = 0; j < Mm; j++){
            Celula[i][j].contador = 0;
            Celula[i][j].lista = new int[N_VIZINHOS_MAX];
            for(int k = 0; k < N_VIZINHOS_MAX; k++)
                Celula[i][j].lista[k] = 0;
        }

    for (int i = 0; i < numero_graos; i++){
        Graos[i].vizinhos[0] = 0;
        int x = (int) ((Graos[i].Linear.Posicao.x -minx) / l);
        int y = (int) ((Graos[i].Linear.Posicao.y -miny) / m);
        x = (x < 0) ? 0 : x;
		x = (x >= Ll) ? Ll -1 : x;
		y = (y < 0) ? 0 : y;
		y = (y >= Mm) ? Mm -1 : y;
        Celula[x][y].lista[Celula[x][y].contador++] = i;
    }
// Neste ponto os graos foram colocados em cada uma das celulas do procura de vizinhos 2

// Agora, realiza-se a busca apenas com os vizinhos da propria caixa e com as caixas vizinhas
    for (int x = 0; x < Ll; x++){
		for (int y = 0; y < Mm; y++){
			for (int k = -1; (k < 2) && (x < Ll); k++){
//                if((x + k >= 0) && (x + k < Ll)){
                {
                    for (int h = 0; (h < 2) && (y < Mm); h++){
//                        if((y+h >= 0) && (y+h < Mm) && !((h== 0) && (k == -1))){
                        if (!((k == -1) && (h == 0))){
                            int d, e;
                            d = x+k;
                            d = d < 1 ? Ll -2 : d; // Condition to the left side, taking the last x position in cell
                            d = d >= Ll-1 ? 1 : d; // Condition to the right side, taking the first x position in cell
                            e = y+h;
                            e = e < 1 ? Mm -2 : e; // Condition to the bottom side, taking the higher y position in cell, as it is builded, it should never happen
                            e = e >= Mm-1 ? 1 : e; // Condition to the top side, taking the lower y position in cell
                            for (int g = 0; g < Celula[x][y].contador; g++){
                                int i = Celula[x][y].lista[g];
                                for (int f = ((h == 0) && (k == 0)) ? g+1 : 0; f < Celula[d][e].contador; f++){
                                    int j = Celula[d][e].lista[f];
                                    double dist = 0.0;
                                    if(CONTROLE == 0){
                                        dist = Graos[i].Distancia(Graos[j]) - Graos[i].Geometria.raio -Graos[j].Geometria.raio;
                                    } else {
                                        dist = Condicao_Periodica_x(Graos[j].Linear.Posicao.x-Graos[j].Linear.Posicao.x);
                                        dist = sqrt(dist*dist +(Graos[j].Linear.Posicao.y-Graos[j].Linear.Posicao.y)*(Graos[j].Linear.Posicao.y-Graos[j].Linear.Posicao.y));
                                        dist -= Graos[i].Geometria.raio +Graos[j].Geometria.raio;
                                    }
                                    if (dist < RAIO_ATUACAO * raio_maior){
                                        Graos[i].vizinhos[++Graos[i].vizinhos[0]] = j;
                                    }
                                }
                            }
                        }
                    }
                }
			}
		}
	}
// Definiu-se os vizinhos graos-graos
	int x = (int) ((Intruso.Linear.Posicao.x -minx) / l);
    int y = (int) ((Intruso.Linear.Posicao.y -miny) / m);
/*    x = (x < 0) ? 0 : x;
    x = (x >= Ll) ? Ll -1 : x;
    y = (y < 0) ? 0 : y;
    y = (y >= Mm) ? Mm -1 : y;

	for (int k = -1; (k < 2) && (k < Ll); k++){
		if ((x+k >= 0) && (x+k < Ll)){
    		for (int h = -1; (h < 2) && (h < Mm); h++){
				if ((y+h >= 0) && (y+h < Mm)){
					for (int g =0; g < Celula[x+k][y+h].contador; g++){
						int i = Celula[x+k][y+h].lista[g];
                        if (Graos[i].Distancia(Intruso) - Graos[i].Geometria.raio -Intruso.Geometria.raio < RAIO_ATUACAO * Intruso.Geometria.raio){
                            Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -1;
                        }
					}
				}
			}
		}
	}

    if (x <= 1) Intruso.vizinhos[++Intruso.vizinhos[0]] = -2;
    if (x >= Ll -2) Intruso.vizinhos[++Intruso.vizinhos[0]] = -3;
    if (y <= 1) Intruso.vizinhos[++Intruso.vizinhos[0]] = -4;
    if (y >= Mm) Intruso.vizinhos[++Intruso.vizinhos[0]] = -5;*/

    for (int i = 0; i < numero_graos; i++){
        double dist = 0.0;
        if(CONTROLE == 0){
            dist = Graos[i].Distancia(Intruso) - Graos[i].Geometria.raio -Intruso.Geometria.raio;
        } else {
            dist = Condicao_Periodica_x(Intruso.Linear.Posicao.x-Graos[i].Linear.Posicao.x);
            dist = sqrt(dist*dist +(Intruso.Linear.Posicao.y-Graos[i].Linear.Posicao.y)*(Intruso.Linear.Posicao.y-Graos[i].Linear.Posicao.y));
            dist -= Intruso.Geometria.raio +Graos[i].Geometria.raio;
        }
        if (dist < RAIO_ATUACAO * Intruso.Geometria.raio){
            Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -1;
        }
	}
	if (Intruso.Linear.Posicao.x - Bordas[0].Linear.Posicao.x -Intruso.Geometria.raio -Bordas[0].Geometria.raio < 1.5 * RAIO_ATUACAO * raio_maior){
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            Intruso.vizinhos[++Intruso.vizinhos[0]] = -2;
        }
    }
    if (Bordas[3].Linear.Posicao.x - Intruso.Linear.Posicao.x -Intruso.Geometria.raio -Bordas[3].Geometria.raio < 1.5 * RAIO_ATUACAO * raio_maior){
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            Intruso.vizinhos[++Intruso.vizinhos[0]] = -3;
        }
    }
    if(CONTROLE==0){
        if (Intruso.Linear.Posicao.y - Bordas[0].Linear.Posicao.y -Intruso.Geometria.raio -Bordas[0].Geometria.raio < 1.5 * RAIO_ATUACAO * raio_maior){
    #ifdef _OPENMP
            #pragma omp critical
    #endif
            {
                Intruso.vizinhos[++Intruso.vizinhos[0]] = -4;
            }
        }
        if (Bordas[3].Linear.Posicao.y - Intruso.Linear.Posicao.y -Intruso.Geometria.raio -Bordas[3].Geometria.raio < 1.5 * RAIO_ATUACAO * raio_maior){
    #ifdef _OPENMP
            #pragma omp critical
    #endif
            {
                Intruso.vizinhos[++Intruso.vizinhos[0]] = -5;
            }
        }
    }

// Definiu-se os vizinhos graos-intruso
    if(CONTROLE==0){
        for (x = 0; x < Ll; x++){
		    for (y = 0; (y < 2) && (y < Mm); y++){
			    for (int g = 0; g < Celula[x][y].contador; g++){
				    int i = Celula[x][y].lista[g];
				    if (Graos[i].Linear.Posicao.y - Bordas[0].Linear.Posicao.y -Graos[i].Geometria.raio -Bordas[0].Geometria.raio < RAIO_ATUACAO * raio_maior){
                        Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -4;
                    }
                }
			    for (int g = 0; g < Celula[x][Mm-1-y].contador; g++){
				    int i = Celula[x][Mm-1-y].lista[g];
				    if (Bordas[3].Linear.Posicao.y - Graos[i].Linear.Posicao.y -Graos[i].Geometria.raio -Bordas[3].Geometria.raio < RAIO_ATUACAO * raio_maior){
                        Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -5;
				    }
			    }
		    }
	    }
	}
// Definiu-se os vizinhos grao-parede X
	for (y = 0; y < Mm; y++){
		for (x  = 0; (x < 2) && (x < Ll); x++){
			for (int g = 0; g < Celula[x][y].contador; g++){
				int i = Celula[x][y].lista[g];
				if (Graos[i].Linear.Posicao.x - Bordas[0].Linear.Posicao.x -Graos[i].Geometria.raio -Bordas[0].Geometria.raio < RAIO_ATUACAO * raio_maior){
                    Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -2;
				}
			}
			for (int g = 0; g < Celula[Ll-1-x][y].contador; g++){
				int i = Celula[Ll-1-x][y].lista[g];
				if (Bordas[3].Linear.Posicao.x - Graos[i].Linear.Posicao.x -Graos[i].Geometria.raio -Bordas[3].Geometria.raio < RAIO_ATUACAO * raio_maior){
                    Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -3;
				}
			}
		}
	}
// Definiu-se os vizinhos grao-parede Y

    for (int i = 0; i < Ll; i++)
        for (int j = 0; j < Mm; j++)
            if (Celula[i][j].lista != NULL)
                delete [] Celula[i][j].lista;

    return 0;
}

int PREDITOR_CORRETOR::Procura_Vizinhos(unsigned int CONTROLE){
    Intruso.vizinhos[0] = 0;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numero_graos; i++){
        Graos[i].vizinhos[0] = 0;
        double dist = 0.0;
        for (int j = i+1; j < numero_graos; j++){
            if(CONTROLE == 0){
                dist = Graos[i].Distancia(Graos[j]) - Graos[i].Geometria.raio -Graos[j].Geometria.raio;
            } else {
                dist = Condicao_Periodica_x(Graos[j].Linear.Posicao.x-Graos[j].Linear.Posicao.x);
                dist = sqrt(dist*dist +(Graos[j].Linear.Posicao.y-Graos[j].Linear.Posicao.y)*(Graos[j].Linear.Posicao.y-Graos[j].Linear.Posicao.y));
                dist -= Graos[i].Geometria.raio -Graos[j].Geometria.raio;
            }
            if (dist < RAIO_ATUACAO * Intruso.Geometria.raio){
                Graos[i].vizinhos[++Graos[i].vizinhos[0]] = j;
            }
        }
        if(CONTROLE == 0){
            dist = Graos[i].Distancia(Intruso) - Graos[i].Geometria.raio -Intruso.Geometria.raio;
        } else {
            dist = Condicao_Periodica_x(Intruso.Linear.Posicao.x-Graos[i].Linear.Posicao.x);
            dist = sqrt(dist*dist +(Intruso.Linear.Posicao.y-Graos[i].Linear.Posicao.y)*(Intruso.Linear.Posicao.y-Graos[i].Linear.Posicao.y));
            dist -= Intruso.Geometria.raio -Graos[i].Geometria.raio;
        }
        if (dist < RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -1;
            }
        }
        if (Graos[i].Linear.Posicao.x - Bordas[0].Linear.Posicao.x -Graos[i].Geometria.raio -Bordas[0].Geometria.raio < RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -2;
            }
        }
        if (Bordas[3].Linear.Posicao.x - Graos[i].Linear.Posicao.x -Graos[i].Geometria.raio -Bordas[3].Geometria.raio < RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -3;
            }
        }
        if (CONTROLE == 0){
            if (Graos[i].Linear.Posicao.y - Bordas[0].Linear.Posicao.y -Graos[i].Geometria.raio -Bordas[0].Geometria.raio < RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -4;
                }
            }
            if (Bordas[3].Linear.Posicao.y - Graos[i].Linear.Posicao.y -Graos[i].Geometria.raio -Bordas[0].Geometria.raio < RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    Graos[i].vizinhos[++Graos[i].vizinhos[0]] = -5;
                }
            }
        }
    }

    if (Intruso.Linear.Posicao.x - Bordas[0].Linear.Posicao.x -Intruso.Geometria.raio -Bordas[0].Geometria.raio < 1.5 * RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            Intruso.vizinhos[++Intruso.vizinhos[0]] = -2;
        }
    }
    if (Bordas[3].Linear.Posicao.x - Intruso.Linear.Posicao.x -Intruso.Geometria.raio -Bordas[3].Geometria.raio < 1.5 * RAIO_ATUACAO * Intruso.Geometria.raio){
#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            Intruso.vizinhos[++Intruso.vizinhos[0]] = -3;
        }
    }
    if (CONTROLE == 0){
        if (Intruso.Linear.Posicao.y - Bordas[0].Linear.Posicao.y -Intruso.Geometria.raio -Bordas[0].Geometria.raio < 1.5 * RAIO_ATUACAO * Intruso.Geometria.raio){
    #ifdef _OPENMP
            #pragma omp critical
    #endif
            {
                Intruso.vizinhos[++Intruso.vizinhos[0]] = -4;
            }
        }
        if (Bordas[3].Linear.Posicao.y - Intruso.Linear.Posicao.y -Intruso.Geometria.raio -Bordas[3].Geometria.raio < 1.5 * RAIO_ATUACAO * Intruso.Geometria.raio){
    #ifdef _OPENMP
            #pragma omp critical
    #endif
            {
                Intruso.vizinhos[++Intruso.vizinhos[0]] = -5;
            }
        }
    }

    return 0;
}

int PREDITOR_CORRETOR::Preditor(unsigned int parametro){
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numero_graos; i++){
        if (parametro / 2 == 0){
            Graos[i].Linear.Posicao.x += Graos[i].Linear.Velocidade.x*dt + Graos[i].Linear.Aceleracao.x*dt*dt/2;
        } else {
            Graos[i].Linear.Posicao.x = Condicao_Periodica_x(Graos[i].Linear.Posicao.x +Graos[i].Linear.Velocidade.x*dt + Graos[i].Linear.Aceleracao.x*dt*dt/2);
        }
        Graos[i].Linear.Posicao.y += Graos[i].Linear.Velocidade.y*dt + Graos[i].Linear.Aceleracao.y*dt*dt/2;
        Graos[i].Angular.Posicao.x += Graos[i].Angular.Velocidade.x*dt + Graos[i].Angular.Aceleracao.x*dt*dt/2;

        Graos[i].PrevistoL.Aceleracao = {Graos[i].Linear.Aceleracao.x, Graos[i].Linear.Aceleracao.y};
        Graos[i].PrevistoA.Aceleracao.x = Graos[i].Angular.Aceleracao.x; // x nesse caso representa alpha
        Graos[i].PrevistoL.Velocidade = {Graos[i].Linear.Velocidade.x +Graos[i].PrevistoL.Aceleracao.x*dt, Graos[i].Linear.Velocidade.y +Graos[i].PrevistoL.Aceleracao.y*dt};
        Graos[i].PrevistoA.Velocidade.x = Graos[i].Angular.Velocidade.x +Graos[i].PrevistoA.Aceleracao.x*dt; //x nesse caso representa omega
    }
    if (parametro % 2){
        if (parametro / 2 == 0){
            Intruso.Linear.Posicao.x += Intruso.Linear.Velocidade.x*dt + Intruso.Linear.Aceleracao.x*dt*dt/2;
        } else {
            Intruso.Linear.Posicao.x = Condicao_Periodica_x(Intruso.Linear.Posicao.x +Intruso.Linear.Velocidade.x*dt + Intruso.Linear.Aceleracao.x*dt*dt/2);
        }
        Intruso.Linear.Posicao.y += Intruso.Linear.Velocidade.y*dt + Intruso.Linear.Aceleracao.y*dt*dt/2;
        Intruso.Angular.Posicao.x += Intruso.Angular.Velocidade.x*dt + Intruso.Angular.Aceleracao.x*dt*dt/2;

        Intruso.PrevistoL.Aceleracao = {Intruso.Linear.Aceleracao.x, Intruso.Linear.Aceleracao.y};
        Intruso.PrevistoA.Aceleracao.x = Intruso.Angular.Aceleracao.x; // x nesse caso representa alpha
        Intruso.PrevistoL.Velocidade = {Intruso.Linear.Velocidade.x +Intruso.PrevistoL.Aceleracao.x*dt, Intruso.Linear.Velocidade.y +Intruso.PrevistoL.Aceleracao.y*dt};
        Intruso.PrevistoA.Velocidade.x = Intruso.Angular.Velocidade.x +Intruso.PrevistoA.Aceleracao.x*dt; //x nesse caso representa omega
    }
    return 0;
}

int PREDITOR_CORRETOR::Subir_Parede(double v){
    for(int i = 0; i < 4; i++)
        Bordas[i].Linear.Posicao.y += v*dt;

    return 0;
}

int PREDITOR_CORRETOR::Vibra_Parede(double A, double T){
    for(int i = 0; i < 4; i++)
        Bordas[i].Linear.Posicao.y += A * sin(2*PI*passos/T)*(2*PI/T);

    return 0;
}

int PREDITOR_CORRETOR::Detectar_Contatos(unsigned int CONTROLE){
    Intruso.contatos[0] = 0;

    //Recuperacao das forcas tangenciais
    REACAO_TANGENCIAL RT[numero_graos +1]; // Reacao tangencial do passo de tempo anterior
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i <= numero_graos; i++){
        RT[i].forca = NULL;
        RT[i].contato = NULL;
        RT[i].forca = new double[N_CONTATOS_MAX];
        RT[i].contato = new int[N_CONTATOS_MAX];
        if ((RT[i].forca == NULL) || (RT[i].contato == NULL))
            fprintf(stderr, "Erro ao alocar memoria na RT!\n");
    }
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numero_graos; i++){
        RT[i].contato[0] = Graos[i].Reacao_Tangencial.contato[0];
        RT[i].forca[0] = 0.0;
        for (int j = 1; j <= Graos[i].Reacao_Tangencial.contato[0]; j++){
            RT[i].forca[j] = Graos[i].Reacao_Tangencial.forca[j];
            RT[i].contato[j] = Graos[i].Reacao_Tangencial.contato[j];
            Graos[i].Reacao_Tangencial.forca[j] = 0;
            Graos[i].Reacao_Tangencial.contato[j] = -10;
        }
        for (int j = Graos[i].Reacao_Tangencial.contato[0] +1; j < N_CONTATOS_MAX; j++){
            RT[i].forca[j] = 0.0;
            RT[i].contato[j] = -10;
            Graos[i].Reacao_Tangencial.forca[j] = 0;
            Graos[i].Reacao_Tangencial.contato[j] = -10;
        }
    }
    RT[numero_graos].contato[0] = Intruso.Reacao_Tangencial.contato[0];
    RT[numero_graos].forca[0] = 0;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int j = 1; j <= Intruso.Reacao_Tangencial.contato[0]; j++){
        RT[numero_graos].forca[j] = Intruso.Reacao_Tangencial.forca[j];
        RT[numero_graos].contato[j] = Intruso.Reacao_Tangencial.contato[j];
        Intruso.Reacao_Tangencial.forca[j] = 0;
        Intruso.Reacao_Tangencial.contato[j] = -10;
    }
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int j = Intruso.Reacao_Tangencial.contato[0] +1; j < N_CONTATOS_MAX; j++){
        RT[numero_graos].forca[j] = 0.0;
        RT[numero_graos].contato[j] = -10;
        Intruso.Reacao_Tangencial.forca[j] = 0;
        Intruso.Reacao_Tangencial.contato[j] = -10;
    }
// Zerou as Reacoes Tangenciais no passo atual
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = 0; i < numero_graos; i++){
        double interpenetracao = 0.0;
        Graos[i].contatos[0] = 0;
        for(int j = 1; j <= Graos[i].vizinhos[0]; j++){
			if (Graos[i].vizinhos[j] >= 0){
			    if (CONTROLE == 0){
                    interpenetracao = Graos[i].Interpenetracao(Graos[Graos[i].vizinhos[j]]);
                } else {
                    interpenetracao = sqrt(pow(Condicao_Periodica_x(Graos[Graos[i].vizinhos[j]].Linear.Posicao.x -Graos[i].Linear.Posicao.x), 2) +pow(Graos[Graos[i].vizinhos[j]].Linear.Posicao.y -Graos[i].Linear.Posicao.y, 2));
                }
            }
            if (Graos[i].vizinhos[j] == -1){
                if (CONTROLE == 0){
                    interpenetracao = Graos[i].Interpenetracao(Intruso);
                } else {
                    interpenetracao = sqrt(pow(Condicao_Periodica_x(Intruso.Linear.Posicao.x -Graos[i].Linear.Posicao.x), 2) +pow(Intruso.Linear.Posicao.y -Graos[i].Linear.Posicao.y, 2));
                }
            }
            if (Graos[i].vizinhos[j] == -2){
                interpenetracao = ((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.x - Graos[i].Geometria.raio));
            }
            if (Graos[i].vizinhos[j] == -3){
                interpenetracao = -((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.x + Graos[i].Geometria.raio));
            }
            if (CONTROLE == 0){
                if (Graos[i].vizinhos[j] == -4){
                    interpenetracao = ((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.y - Graos[i].Geometria.raio));
                }
                if (Graos[i].vizinhos[j] == -5){
                    interpenetracao = -((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.y + Graos[i].Geometria.raio));
                }
            }
            if (interpenetracao > 0.0){
#ifdef _OPENMP
                #pragma omp atomic
#endif
                Graos[i].contatos[0]++;
                Graos[i].contatos[Graos[i].contatos[0]] = Graos[i].vizinhos[j];
                Graos[i].Reacao_Tangencial.contato[Graos[i].contatos[0]] = Graos[i].vizinhos[j];
                int k = 0;
                for (k = 1; (RT[i].contato[k] != Graos[i].vizinhos[j]) && (RT[i].contato[k] != -10) && (k < N_CONTATOS_MAX); k++); // Busca da forca tangencial do passo anterior
                Graos[i].Reacao_Tangencial.forca[Graos[i].contatos[0]] = (RT[i].contato[k] != -10) ? RT[i].forca[k] : 0; // Atualizando a lista de reacoes tangenciais
            }
        }
        Graos[i].Reacao_Tangencial.contato[0] = Graos[i].contatos[0];
    }
// Encontrou o contato dos graos

// Encontrando as bordas para a condicao de contato com o Intruso
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int j = 1; j <= Intruso.vizinhos[0]; j++){
        double interpenetracao = 0.0;

        if (Intruso.vizinhos[j] == -2){
            interpenetracao = ((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.x - Intruso.Geometria.raio));
        }
        if (Intruso.vizinhos[j] == -3){
            interpenetracao = -((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.x + Intruso.Geometria.raio));
        }
        if (CONTROLE == 0){
            if (Intruso.vizinhos[j] == -4){
                interpenetracao = ((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.y - Intruso.Geometria.raio));
            }
            if (Intruso.vizinhos[j] == -5){
                interpenetracao = -((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.y + Intruso.Geometria.raio));
            }
        }
        if (interpenetracao > 0.0){
#ifdef _OPENMP
            #pragma omp atomic
#endif
            Intruso.contatos[0]++;
            Intruso.contatos[Intruso.contatos[0]] = Intruso.vizinhos[j];
            Intruso.Reacao_Tangencial.contato[Intruso.contatos[0]] = Intruso.vizinhos[j];
            int k = 0;
            for (k = 1; (RT[numero_graos].contato[k] != Intruso.vizinhos[j]) && (RT[numero_graos].contato[k] != -10) && (k < N_CONTATOS_MAX); k++); // Busca da forca tangencial do passo anterior
            Intruso.Reacao_Tangencial.forca[Intruso.contatos[0]] = (RT[numero_graos].contato[k] != -10) ? RT[numero_graos].forca[k] : 0; // Atualizando a lista de reacoes tangenciais
        }
    }
    Intruso.Reacao_Tangencial.contato[0] = Intruso.contatos[0];
// Encontrou os contatos do Intruso
    for (int i = 0; i <= numero_graos; i++){
        if (RT[i].forca != NULL){
            delete [] RT[i].forca;
            RT[i].forca = NULL;
        }
        if (RT[i].contato != NULL){
            delete [] RT[i].contato;
            RT[i].contato = NULL;
        }
    }

    return 0;
}

int PREDITOR_CORRETOR::Calculo_Forca(unsigned int CONTROLE){
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numero_graos; i++){
        Graos[i].Forca.x = 0.0;
        Graos[i].Forca.y = -gravidade *Graos[i].massa;
        Graos[i].Torque = 0.0;
    }
    Intruso.Forca.x = 0.0;
    Intruso.Forca.y = -gravidade *Intruso.massa;
    Intruso.Torque = 0.0;

    for (int i = 0; i < 4; i++){
        Bordas[i].Forca.x = 0.0;
        Bordas[i].Forca.y = 0.0;
        Bordas[i].Torque = 0.0;
    }

#ifdef _OPENMP
    #pragma omp parallel for ordered
#endif
    for (int i = 0; i < numero_graos; i++){
        VETOR Velocidade[2], RI; //RI -> Razao de interpenetracao
        double forca_normal, forca_tangencial, modulo, interpenetracao, atrito;
        GRAO *G1 = NULL, *G2 = NULL;

        G1 = &Graos[i];
        Velocidade[0] = G1->PrevistoL.Velocidade;
        Velocidade[1] = {0.0, 0.0};

        for (int j = 1; j <= G1->contatos[0]; j++){
            if (Graos[i].contatos[j] >= 0){
                G2 = &Graos[Graos[i].contatos[j]];
                Velocidade[1] = G2->PrevistoL.Velocidade;

                if (CONTROLE == 0){
                    RI.x = G2->Linear.Posicao.x -G1->Linear.Posicao.x;
                } else {
                    RI.x = Condicao_Periodica_x(G2->Linear.Posicao.x -G1->Linear.Posicao.x);
                }
                RI.y = G2->Linear.Posicao.y -G1->Linear.Posicao.y;
                modulo = sqrt(RI.x*RI.x + RI.y*RI.y);
                RI.x /= modulo;
                RI.y /= modulo;

                if (CONTROLE == 0){
                    interpenetracao = G1->Interpenetracao(*G2);
                } else {
                    interpenetracao = sqrt(pow(Condicao_Periodica_x(G1->Linear.Posicao.x -G2->Linear.Posicao.x), 2) +pow(G1->Linear.Posicao.y -G2->Linear.Posicao.y, 2));
                }
                atrito = G1->atrito[0] + G2->atrito[0];
            }
            if (G1->contatos[j] == -1){
                G2 = &Intruso;
                Velocidade[1] = G2->PrevistoL.Velocidade;

                if (CONTROLE == 0){
                    RI.x = G2->Linear.Posicao.x -G1->Linear.Posicao.x;
                } else {
                    RI.x = Condicao_Periodica_x(G2->Linear.Posicao.x -G1->Linear.Posicao.x);
                }
                RI.y = G2->Linear.Posicao.y - G1->Linear.Posicao.y;
                modulo = sqrt(RI.x*RI.x +  RI.y*RI.y);
                RI.x /= modulo;
                RI.y /= modulo;

                if (CONTROLE == 0){
                    interpenetracao = G1->Interpenetracao(*G2);
                } else {
                    interpenetracao = sqrt(pow(Condicao_Periodica_x(G1->Linear.Posicao.x -G2->Linear.Posicao.x), 2) +pow(G1->Linear.Posicao.y -G2->Linear.Posicao.y, 2));
                }
                atrito = G1->atrito[1] + G2->atrito[1];
            }
            if (G1->contatos[j] == -2){
                G2 = &Bordas[0];

                RI.x = -fabs(G2->Linear.Posicao.x -G1->Linear.Posicao.x);
                RI.y = 0;
                modulo = fabs(RI.x);
                RI.x /= modulo;

                interpenetracao = ((G2->Linear.Posicao.x + G2->Geometria.raio) - (G1->Linear.Posicao.x - G1->Geometria.raio));
                atrito = G1->atrito[1] + G2->atrito[1];
            }
            if (G1->contatos[j] == -3){
                G2 = &Bordas[3];

                RI.x = +fabs(G2->Linear.Posicao.x -G1->Linear.Posicao.x);
                RI.y = 0;
                modulo = fabs(RI.x);
                RI.x /= modulo;

                interpenetracao = -((G2->Linear.Posicao.x - G2->Geometria.raio) - (G1->Linear.Posicao.x + G1->Geometria.raio));
                atrito = G1->atrito[1] + G2->atrito[1];
            }
            if (CONTROLE == 0){
                if (G1->contatos[j] == -4){
                    G2 = &Bordas[0];

                    RI.x = 0;
                    RI.y = -fabs(G2->Linear.Posicao.y -G1->Linear.Posicao.y);
                    modulo = fabs(RI.y);
                    RI.y /= modulo;

                    interpenetracao = ((G2->Linear.Posicao.y + G2->Geometria.raio) - (G1->Linear.Posicao.y - G1->Geometria.raio));
                    atrito = G1->atrito[1] + G2->atrito[1];
                }
                if (G1->contatos[j] == -5){
                    G2 = &Bordas[3];

                    RI.x = 0;
                    RI.y = +fabs(G2->Linear.Posicao.y -G1->Linear.Posicao.y);
                    modulo = fabs(RI.y);
                    RI.y /= modulo;

                    interpenetracao = -((G2->Linear.Posicao.y - G2->Geometria.raio) - (G1->Linear.Posicao.y + G1->Geometria.raio));
                    atrito = G1->atrito[1] + G2->atrito[1];
                }
            }
//CALCULO EFETIVO DAS FORCAS
            double kn = (G1->mola.normal + G2->mola.normal)/2; // Mola na direcao normal associada em paralelo (Deformacao fixa!)
            double kt = (G1->mola.tangencial + G2->mola.tangencial)/2;
            forca_normal = interpenetracao * kn; // Parcela Elastica
            forca_normal -= (RI.x * (Velocidade[1].x -Velocidade[0].x) + RI.y * (Velocidade[1].y -Velocidade[0].y)) * gn; // Parcela Viscosa

            forca_tangencial = G1->Reacao_Tangencial.forca[j];
            forca_tangencial -= (-RI.y * (Velocidade[1].x - Velocidade[0].x) + RI.x * (Velocidade[1].y - Velocidade[0].y) - G1->PrevistoA.Velocidade.x * G1->Geometria.raio - G2->PrevistoA.Velocidade.x * G2->Geometria.raio) * kt * dt;

            atrito = 0.5;
            double atrito_maximo = interpenetracao * kn * atrito;
            if (fabs(forca_tangencial) > atrito_maximo){
                forca_tangencial = (forca_tangencial > 0) ? atrito_maximo : -atrito_maximo;
            }

#ifdef _OPENMP
            #pragma omp ordered
#endif
            {
                G1->Forca.x -= forca_normal * RI.x - forca_tangencial * RI.y;
                G1->Forca.y -= forca_normal * RI.y + forca_tangencial * RI.x;

                G2->Forca.x += forca_normal * RI.x - forca_tangencial * RI.y;
                G2->Forca.y += forca_normal * RI.y + forca_tangencial * RI.x;

                G1->Torque -= forca_tangencial * G1->Geometria.raio;
                G2->Torque -= forca_tangencial * G2->Geometria.raio;

                G1->Reacao_Tangencial.forca[j] = forca_tangencial;
            }
        }
    }

    // Forca do Intruso com as paredes
#ifdef _OPENMP
    #pragma omp parallel for ordered
#endif
    for (int j = 1; j <= Intruso.contatos[0]; j++){
        VETOR Velocidade[1], RI; //RI -> Razao de interpenetracao
        double forca_normal, forca_tangencial, modulo, interpenetracao, atrito;
        GRAO *G1 = NULL, *G2 = NULL;

        G1 = &Intruso;
        Velocidade[0] = G1->PrevistoL.Velocidade;
        Velocidade[1] = {0.0, 0.0};

        if (G1->contatos[j] == -2){
            G2 = &Bordas[0];

            RI.x = -fabs(G2->Linear.Posicao.x -G1->Linear.Posicao.x);
            RI.y = 0;
            modulo = fabs(RI.x);
            RI.x /= modulo;

            interpenetracao = fabs((G2->Geometria.raio + G1->Geometria.raio + G2->Linear.Posicao.x - G1->Linear.Posicao.x ));
            atrito = G1->atrito[2] + G2->atrito[2];
        }
        if (G1->contatos[j] == -3){
            G2 = &Bordas[3];

            RI.x = +fabs(G2->Linear.Posicao.x -G1->Linear.Posicao.x);
            RI.y = 0;
            modulo = fabs(RI.x);
            RI.x /= modulo;

            interpenetracao = fabs((G2->Geometria.raio + G1->Geometria.raio - G2->Linear.Posicao.x + G1->Linear.Posicao.x ));
            atrito = G1->atrito[2] + G2->atrito[2];
        }
        if (CONTROLE == 0){
            if (G1->contatos[j] == -4){
                G2 = &Bordas[0];

                RI.x = 0;
                RI.y = -fabs(G2->Linear.Posicao.y -G1->Linear.Posicao.y);
                modulo = fabs(RI.y);
                RI.y /= modulo;

                interpenetracao = fabs((G2->Geometria.raio + G1->Geometria.raio + G2->Linear.Posicao.y - G1->Linear.Posicao.y ));
                atrito = G1->atrito[2] + G2->atrito[2];
            }
            if (G1->contatos[j] == -5){
                G2 = &Bordas[3];

                RI.x = 0;
                RI.y = +fabs(G2->Linear.Posicao.y -G1->Linear.Posicao.y);
                modulo = fabs(RI.y);
                RI.y /= modulo;

                interpenetracao = fabs((G2->Geometria.raio + G1->Geometria.raio - G2->Linear.Posicao.y + G1->Linear.Posicao.y ));
                atrito = G1->atrito[2] + G2->atrito[2];
            }
        }
//CALCULO EFETIVO DAS FORCAS
        double kn = (G1->mola.normal + G2->mola.normal)/2;
        double kt = (G1->mola.tangencial + G2->mola.tangencial)/2;
        forca_normal = interpenetracao * kn; // Parcela Elastica
        forca_normal -= (RI.x * (Velocidade[1].x -Velocidade[0].x) + RI.y * (Velocidade[1].y -Velocidade[0].y)) * gn; // Parcela Viscosa

        forca_tangencial = G1->Reacao_Tangencial.forca[j];
        forca_tangencial -= (-RI.y * (Velocidade[1].x - Velocidade[0].x) + RI.x * (Velocidade[1].y - Velocidade[0].y) - G1->PrevistoA.Velocidade.x * G1->Geometria.raio - G2->PrevistoA.Velocidade.x * G2->Geometria.raio) * kt * dt;

        atrito = 0.5;
        double atrito_maximo = interpenetracao * kn * atrito;
        if (fabs(forca_tangencial) > atrito_maximo){
            forca_tangencial = (forca_tangencial > 0) ? atrito_maximo : -atrito_maximo;
        }

#ifdef _OPENMP
        #pragma omp ordered
#endif
        {
            G1->Forca.x -= forca_normal * RI.x - forca_tangencial * RI.y;
            G1->Forca.y -= forca_normal * RI.y + forca_tangencial * RI.x;

            G2->Forca.x += forca_normal * RI.x - forca_tangencial * RI.y;
            G2->Forca.y += forca_normal * RI.y + forca_tangencial * RI.x;

            G1->Torque -= forca_tangencial * G1->Geometria.raio;
            G2->Torque -= forca_tangencial * G2->Geometria.raio;

            G1->Reacao_Tangencial.forca[j] = forca_tangencial;
        }
    }

    for (int i = 0; i < 4; i++){
        Bordas[i].Forca.x = 0.0;
        Bordas[i].Forca.y = 0.0;
        Bordas[i].Torque = 0.0;
    }

    return 0;
}

int PREDITOR_CORRETOR::Corretor(unsigned int parametro){
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numero_graos; i++){
        Graos[i].Linear.Aceleracao.x = Graos[i].Forca.x / Graos[i].massa;
        Graos[i].Linear.Aceleracao.y = Graos[i].Forca.y / Graos[i].massa;
        Graos[i].Linear.Velocidade.x = Graos[i].PrevistoL.Velocidade.x + dt/2 * (Graos[i].Linear.Aceleracao.x -Graos[i].PrevistoL.Aceleracao.x);
        Graos[i].Linear.Velocidade.y = Graos[i].PrevistoL.Velocidade.y + dt/2 * (Graos[i].Linear.Aceleracao.y -Graos[i].PrevistoL.Aceleracao.y);
        Graos[i].Angular.Aceleracao.x = Graos[i].Torque / Graos[i].inercia;
        Graos[i].Angular.Velocidade.x = Graos[i].PrevistoA.Velocidade.x + dt/2 * (Graos[i].Angular.Aceleracao.x -Graos[i].PrevistoA.Aceleracao.x);
    }
	if (parametro % 2 == 1){
        Intruso.Linear.Aceleracao.x = Intruso.Forca.x / Intruso.massa;
        Intruso.Linear.Aceleracao.y = Intruso.Forca.y / Intruso.massa;
        Intruso.Linear.Velocidade.x = Intruso.PrevistoL.Velocidade.x + dt/2 * (Intruso.Linear.Aceleracao.x -Intruso.PrevistoL.Aceleracao.x);
        Intruso.Linear.Velocidade.y = Intruso.PrevistoL.Velocidade.y + dt/2 * (Intruso.Linear.Aceleracao.y -Intruso.PrevistoL.Aceleracao.y);
        Intruso.Angular.Aceleracao.x = Intruso.Torque / Intruso.inercia;
        Intruso.Angular.Velocidade.x = Intruso.PrevistoA.Velocidade.x + dt/2 * (Intruso.Angular.Aceleracao.x -Intruso.PrevistoA.Aceleracao.x);
    }

    return 0;
}

double PREDITOR_CORRETOR::Coordenacao(int parametro){
    double coordenacao = 0.0;
    int lista[numero_graos], contador = 0;


    for (int i = 0; i < numero_graos; i++){
        lista[contador] = 0;
        if (parametro == 0){
            if (Graos[i].Linear.Posicao.y < Intruso.Linear.Posicao.y -Intruso.Geometria.raio)
                lista[contador++] = i;
        }
        else
            lista[contador++] = i;
    }
#ifdef _OPENMP
    #pragma omp parallel for reduction(+: coordenacao)
#endif
    for (int i = 0; i < contador; i++){
        for (int j = 1; j <= Graos[lista[i]].contatos[0]; j++)
            coordenacao += (Graos[lista[i]].contatos[j] >= 0) ? 1.0 : 0.0;
    }
    coordenacao += Intruso.contatos[0];
    coordenacao /= (contador+1.0)/2.0;

    return coordenacao;
}

double PREDITOR_CORRETOR::Compactacao(){
    double phi = (Bordas[3].Linear.Posicao.x -Bordas[3].Geometria.raio -Bordas[0].Linear.Posicao.x -Bordas[0].Geometria.raio) * Intruso.Geometria.raio * 4.0, area = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+: area)
#endif
    for (int i = 0; i < numero_graos; i++){
        if ((Graos[i].Linear.Posicao.y -Graos[i].Geometria.raio < Intruso.Linear.Posicao.y +Intruso.Geometria.raio * 2.0) && (Graos[i].Linear.Posicao.y +Graos[i].Geometria.raio > Intruso.Linear.Posicao.y -Intruso.Geometria.raio * 2.0)){
            area += Graos[i].Geometria.raio*Graos[i].Geometria.raio*PI;
        }
    }
    area += Intruso.Geometria.raio*Intruso.Geometria.raio*PI;
    phi = area/phi;

    return phi;
}

double PREDITOR_CORRETOR::Interpenetracao(){
    double area = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+: area)
#endif
    for (int i = 0; i < numero_graos; i++)
        if (Graos[i].Linear.Posicao.y < Intruso.Linear.Posicao.y -Intruso.Geometria.raio)
            area += Area_Interpenetracao(i);

    return area;
}

double PREDITOR_CORRETOR::Pressao_Parede(int b){
    double pressao = 0.0;

#ifdef _OPENMP
    #pragma omp parallel for reduction(+: pressao)
#endif
    for (int i = 0; i < numero_graos; i++){
        for (int j = 1; j < Graos[i].contatos[0]; j++){
            if (Graos[i].contatos[j] == b){
                if (b == -2){
                    pressao += ((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.x - Graos[i].Geometria.raio));
                }
                if (b == -3){
                    pressao += -((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.x + Graos[i].Geometria.raio));
                }
                if (b == -4){
                    pressao += ((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.y - Graos[i].Geometria.raio));
                }
                if (b == -5){
                    pressao += -((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.y + Graos[i].Geometria.raio));
                }
            }
        }
    }
#ifdef _OPENMP
    #pragma omp parallel for reduction(+: pressao)
#endif
    for (int j = 1; j < Intruso.contatos[0]; j++){
        if (Intruso.contatos[j] == b){
            if (b == -2){
                pressao += ((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.x - Intruso.Geometria.raio));
            }
            if (b == -3){
                pressao += ((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.x + Intruso.Geometria.raio));
            }
            if (b == -4){
                pressao += ((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.y - Intruso.Geometria.raio));
            }
            if (b == -5){
                pressao += ((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.y + Intruso.Geometria.raio));
            }
        }
    }

    return pressao;
}

int PREDITOR_CORRETOR::Verifica_Equilibrio(){
    if ((Bordas[3].Linear.Posicao.y - Intruso.Linear.Posicao.y)/(Bordas[3].Linear.Posicao.y - Bordas[0].Linear.Posicao.y) >= 1/4)
        return 1;

    return 0;
}

int PREDITOR_CORRETOR::Salvar_Configuracao(const char *nome){
    FILE *Arquivo;

    Arquivo = fopen(nome, "w");
    long int contatos = 0;

    for (int i = 0; i < numero_graos; i++){
        contatos += Graos[i].contatos[0];
    }
    contatos += Intruso.contatos[0];

    fprintf(Arquivo, "%d %d %d %llu %li\n", numero_graos, numero_graos + 5, 3, passos, contatos);
    for (int i = 0; i < numero_graos; i++){
        fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Graos[i].Linear.Posicao.x, Graos[i].Linear.Posicao.y, Graos[i].Angular.Posicao.x, Graos[i].Linear.Velocidade.x, Graos[i].Linear.Velocidade.y, Graos[i].Angular.Velocidade.x, Graos[i].Linear.Aceleracao.x, Graos[i].Linear.Aceleracao.y, Graos[i].Angular.Aceleracao.x, Graos[i].Geometria.raio);
    }
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Intruso.Linear.Posicao.x, Intruso.Linear.Posicao.y, Intruso.Angular.Posicao.x, Intruso.Linear.Velocidade.x, Intruso.Linear.Velocidade.y, Intruso.Angular.Velocidade.x, Intruso.Linear.Aceleracao.x, Intruso.Linear.Aceleracao.y, Intruso.Angular.Aceleracao.x, Intruso.Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[0].Linear.Posicao.x, Bordas[0].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[0].Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[1].Linear.Posicao.x, Bordas[1].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[1].Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[2].Linear.Posicao.x, Bordas[2].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[2].Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[3].Linear.Posicao.x, Bordas[3].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[3].Geometria.raio);

    for (int i = 0; i < numero_graos; i++){
        for(int j = 1; j <= Graos[i].contatos[0]; j++){
            fprintf(Arquivo, "%i %i %.10e\n", i, Graos[i].contatos[j], Graos[i].Reacao_Tangencial.forca[j]);
        }
    }

    for(int j = 1; j <= Intruso.contatos[0]; j++){
        fprintf(Arquivo, "%i %i %.10e\n", -1, Intruso.contatos[j], Intruso.Reacao_Tangencial.forca[j]);
    }

    fclose(Arquivo);

/*    sprintf(nome, "%s.bin", nome);
    Arquivo = fopen(nome, "wb");

    auto tmp = numero_graos + 5;

    fwrite(&numero_graos, sizeof(numero_graos), 1, Arquivo);
    fwrite(&tmp, sizeof(tmp), 1, Arquivo);
    fwrite(&passos, sizeof(passos), 1, Arquivo);

    // Todo o vetor de grãos em ordem: Raio, Posição, Velocidade, Aceleração

    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Geometria.raio, sizeof(double), 1, Arquivo);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Linear.Posicao.x);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Linear.Posicao.y);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Angular.Posicao.x);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Linear.Velocidade.x);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Linear.Velocidade.y);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Angular.Velocidade.x);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Linear.Aceleracao.x);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Linear.Aceleracao.y);
    }
    for (int i = 0; i < numero_graos; i++){
        fwrite(&Graos[i].Angular.Aceleracao.x);
    }

    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Intruso.Linear.Posicao.x, Intruso.Linear.Posicao.y, Intruso.Angular.Posicao.x, Intruso.Linear.Velocidade.x, Intruso.Linear.Velocidade.y, Intruso.Angular.Velocidade.x, Intruso.Linear.Aceleracao.x, Intruso.Linear.Aceleracao.y, Intruso.Angular.Aceleracao.x, Intruso.Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[0].Linear.Posicao.x, Bordas[0].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[0].Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[1].Linear.Posicao.x, Bordas[1].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[1].Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[2].Linear.Posicao.x, Bordas[2].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[2].Geometria.raio);
    fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Bordas[3].Linear.Posicao.x, Bordas[3].Linear.Posicao.y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Bordas[3].Geometria.raio);

/*	fwrite(&ngtot,sizeof(ngtot),1,out3);
	fwrite(radius,sizeof(double long),ngtot+1,out3);
	fwrite(rx,sizeof(double long),ngtot+1,out3);
	fwrite(ry,sizeof(double long),ngtot+1,out3);
	fwrite(vx,sizeof(double long),ngtot+1,out3);
	fwrite(vy,sizeof(double long),ngtot+1,out3);
	fwrite(ome,sizeof(double long),ngtot+1,out3);
	fwrite(ax,sizeof(double long),ngtot+1,out3);
	fwrite(ay,sizeof(double long),ngtot+1,out3);
	fwrite(domedt,sizeof(double long),ngtot+1,out3);
	fwrite(&ncontfree,sizeof(ncontfree),1,out3);
	fwrite(&ncontot,sizeof(ncontot),1,out3);
	fwrite(ior,sizeof(double long),ncontot+1,out3);
	fwrite(iex,sizeof(double long),ncontot+1,out3);
	fwrite(react,sizeof(double long),ncontot+1,out3);*/

    //fclose(Arquivo);

    return 0;
}

int PREDITOR_CORRETOR::Salvar_Configuracao(){
    return Salvar_Configuracao("termalizado.dat");
}

double PREDITOR_CORRETOR::Energia_Cinetica(){
    double energia = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:energia)
#endif
    for (int i = 0; i < numero_graos; i++){
        energia += Graos[i].massa *(Graos[i].Linear.Velocidade.norma())*(Graos[i].Linear.Velocidade.norma());
    }
    energia += Intruso.massa *(Intruso.Linear.Velocidade.norma())*(Intruso.Linear.Velocidade.norma());
    return energia / 2.0;
}

double PREDITOR_CORRETOR::Energia_Cinetica(int i){
    return Graos[i].massa *(Graos[i].Linear.Velocidade.norma())*(Graos[i].Linear.Velocidade.norma()) / 2.0;
}

double PREDITOR_CORRETOR::Energia_Rotacional(){
    double energia = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+: energia)
#endif
    for (int i = 0; i < numero_graos; i++){
        energia += Graos[i].inercia *(Graos[i].Angular.Velocidade.x)*(Graos[i].Angular.Velocidade.x);
    }
    energia += Intruso.inercia *(Intruso.Angular.Velocidade.x)*(Intruso.Angular.Velocidade.x);
    return energia / 2.0;
}

double PREDITOR_CORRETOR::Energia_Rotacional(int i){
    return Graos[i].inercia *(Graos[i].Angular.Velocidade.x)*(Graos[i].Angular.Velocidade.x);
}

double PREDITOR_CORRETOR::Energia_Elastica(){
    double energia = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction (+: energia)
#endif
    for (int i = 0; i < numero_graos; i++){
        for (int j = 1; j <= Graos[i].contatos[0]; j++){
            if (Graos[i].contatos[j] >= 0)
                energia += (Graos[i].mola.normal+Graos[Graos[i].contatos[j]].mola.normal)/2 *(Graos[i].Interpenetracao(Graos[Graos[i].contatos[j]]))*(Graos[i].Interpenetracao(Graos[Graos[i].contatos[j]]));
            else if (Graos[i].contatos[j] == -1)
                energia += (Graos[i].mola.normal+Intruso.mola.normal)/2 *(Graos[i].Interpenetracao(Intruso))*(Graos[i].Interpenetracao(Intruso));
            else if (Graos[i].contatos[j] == -2)
                energia += (Graos[i].mola.normal+Bordas[0].mola.normal)/2 *2*((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.x - Graos[i].Geometria.raio))*((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.x - Graos[i].Geometria.raio));
            else if (Graos[i].contatos[j] == -3)
                energia += (Graos[i].mola.normal+Bordas[3].mola.normal)/2 *2*((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.x + Graos[i].Geometria.raio))*((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.x + Graos[i].Geometria.raio));
            else if (Graos[i].contatos[j] == -4)
                energia += (Graos[i].mola.normal+Bordas[0].mola.normal)/2 *2*((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.y - Graos[i].Geometria.raio))*((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Graos[i].Linear.Posicao.y - Graos[i].Geometria.raio));
            else if (Graos[i].contatos[j] == -5)
                energia += (Graos[i].mola.normal+Bordas[3].mola.normal)/2 *2*((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.y + Graos[i].Geometria.raio))*((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Graos[i].Linear.Posicao.y + Graos[i].Geometria.raio));
        }
    }
#ifdef _OPENMP
    #pragma omp parallel for reduction (+: energia)
#endif
    for (int j = 1; j <= Intruso.contatos[0]; j++){
        if (Intruso.contatos[j] >= 0)
            energia += (Graos[j].mola.normal+Intruso.mola.normal)/2 *(Graos[j].Interpenetracao(Intruso))*(Graos[j].Interpenetracao(Intruso));
        else if (Intruso.contatos[j] == -2)
            energia += (Intruso.mola.normal+Bordas[0].mola.normal)/2 *2*((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.x - Intruso.Geometria.raio))*((Bordas[0].Linear.Posicao.x + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.x - Intruso.Geometria.raio));
        else if (Intruso.contatos[j] == -3)
            energia += (Intruso.mola.normal+Bordas[3].mola.normal)/2 *2*((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.x + Intruso.Geometria.raio))*((Bordas[3].Linear.Posicao.x - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.x + Intruso.Geometria.raio));
        else if (Intruso.contatos[j] == -4)
            energia += (Intruso.mola.normal+Bordas[0].mola.normal)/2 *2*((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.y - Intruso.Geometria.raio))*((Bordas[0].Linear.Posicao.y + Bordas[0].Geometria.raio) - (Intruso.Linear.Posicao.y - Intruso.Geometria.raio));
        else if (Intruso.contatos[j] == -5)
            energia += (Intruso.mola.normal+Bordas[3].mola.normal)/2 *2*((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.y + Intruso.Geometria.raio))*((Bordas[3].Linear.Posicao.y - Bordas[3].Geometria.raio) - (Intruso.Linear.Posicao.y + Intruso.Geometria.raio));
    }
    return energia / 2.0;
}

double PREDITOR_CORRETOR::Energia_Elastica(int i){
    double energia = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction (+: energia)
#endif
    for (int j = 1; j <= Graos[i].contatos[0]; j++){
        energia += (Graos[i].mola.normal*Graos[Graos[i].contatos[j]].mola.normal)/(Graos[i].mola.normal+Graos[Graos[i].contatos[j]].mola.normal) *(Graos[i].Interpenetracao(Graos[Graos[i].contatos[j]]))*(Graos[i].Interpenetracao(Graos[Graos[i].contatos[j]]));
    }
    return energia / 2.0;
}

double PREDITOR_CORRETOR::Energia_Gravitacional(){
    double energia = 0;
#ifdef _OPENMP
    #pragma omp parallel for reduction (+: energia)
#endif
    for (int i = 0; i < numero_graos; i++){
        energia += Graos[i].massa * (Graos[i].Linear.Posicao.y -Bordas[0].Linear.Posicao.y) * gravidade;
    }
    energia += Intruso.massa * (Intruso.Linear.Posicao.y -Bordas[0].Linear.Posicao.y) * gravidade;
    return energia;
}

double PREDITOR_CORRETOR::Energia_Gravitacional(int i){
    return Graos[i].massa * (Graos[i].Linear.Posicao.y -Bordas[0].Linear.Posicao.y) * gravidade;
}

double PREDITOR_CORRETOR::Area_Interpenetracao(int i){
    double interpenetracao = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for reduction (+: interpenetracao)
#endif
    for (int j = 1; j <= Graos[i].contatos[0]; j++){
        double area = 0;
        if (Graos[i].contatos[j] >= 0){
            double distancia = Graos[i].Distancia(Graos[Graos[i].contatos[j]]);
            double r1 = Graos[i].Geometria.raio;
            double r2 = Graos[Graos[i].contatos[j]].Geometria.raio;
            double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
            double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
            area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
        }
        else if (Graos[i].contatos[j] == -1){
            double distancia = Graos[i].Distancia(Intruso);
            double r1 = Graos[i].Geometria.raio;
            double r2 = Intruso.Geometria.raio;
            double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
            double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
            area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
        }
        else if (Graos[i].contatos[j] == -2){
            double distancia = ((Bordas[0].Linear.Posicao.x) - (Graos[i].Linear.Posicao.x));
            double r1 = Graos[i].Geometria.raio;
            double r2 = Bordas[0].Geometria.raio;
            double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
            double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
            area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
        }
        else if (Graos[i].contatos[j] == -3){
            double distancia = -((Bordas[3].Linear.Posicao.x) - (Graos[i].Linear.Posicao.x));
            double r1 = Graos[i].Geometria.raio;
            double r2 = Bordas[3].Geometria.raio;
            double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
            double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
            area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
        }
        else if (Graos[i].contatos[j] == -4){
            double distancia = ((Bordas[0].Linear.Posicao.y) - (Graos[i].Linear.Posicao.y));
            double r1 = Graos[i].Geometria.raio;
            double r2 = Bordas[3].Geometria.raio;
            double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
            double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
            area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
        }
        else if (Graos[i].contatos[j] == -5){
            double distancia = -((Bordas[3].Linear.Posicao.y) - (Graos[i].Linear.Posicao.y));
            double r1 = Graos[i].Geometria.raio;
            double r2 = Graos[Graos[i].contatos[j]].Geometria.raio;
            double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
            double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
            area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
        }
        interpenetracao += area;
    }
    return interpenetracao;
}

double PREDITOR_CORRETOR::Area_Interpenetracao_Total(){
    double interpenetracao = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for reduction (+: interpenetracao)
#endif
    for (int i = 0; i < numero_graos; i++){
        for (int j = 1; j <= Graos[i].contatos[0]; j++){
            double area = 0;
            if (Graos[i].contatos[j] >= 0){
                double distancia = Graos[i].Distancia(Graos[Graos[i].contatos[j]]);
                double r1 = Graos[i].Geometria.raio;
                double r2 = Graos[Graos[i].contatos[j]].Geometria.raio;
                double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
                double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
                area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
            }
            else if (Graos[i].contatos[j] == -1){
                double distancia = Graos[i].Distancia(Intruso);
                double r1 = Graos[i].Geometria.raio;
                double r2 = Intruso.Geometria.raio;
                double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
                double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
                area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
            }
            else if (Graos[i].contatos[j] == -2){
                double distancia = ((Bordas[0].Linear.Posicao.x) - (Graos[i].Linear.Posicao.x));
                double r1 = Graos[i].Geometria.raio;
                double r2 = Bordas[0].Geometria.raio;
                double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
                double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
                area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
            }
            else if (Graos[i].contatos[j] == -3){
                double distancia = -((Bordas[3].Linear.Posicao.x) - (Graos[i].Linear.Posicao.x));
                double r1 = Graos[i].Geometria.raio;
                double r2 = Bordas[3].Geometria.raio;
                double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
                double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
                area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
            }
            else if (Graos[i].contatos[j] == -4){
                double distancia = ((Bordas[0].Linear.Posicao.y) - (Graos[i].Linear.Posicao.y));
                double r1 = Graos[i].Geometria.raio;
                double r2 = Bordas[3].Geometria.raio;
                double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
                double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
                area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
            }
            else if (Graos[i].contatos[j] == -5){
                double distancia = -((Bordas[3].Linear.Posicao.y) - (Graos[i].Linear.Posicao.y));
                double r1 = Graos[i].Geometria.raio;
                double r2 = Graos[Graos[i].contatos[j]].Geometria.raio;
                double alpha = 2*acos((distancia*distancia+r1*r1-r2*r2)/(2*distancia*r1));
                double beta = 2*acos((distancia*distancia+r2*r2-r1*r1)/(2*distancia*r2));
                area = r1*r1*(alpha -sin(alpha)/2)+r2*r2*(beta -sin(beta)/2);
            }
            interpenetracao += area;
        }
    }
    return interpenetracao;
}
