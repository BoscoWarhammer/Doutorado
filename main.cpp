#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "pushbox.cpp"

#ifdef _OPENMP
    #include <omp.h>
#endif

void help(const char *);

typedef struct Tempo{
    unsigned int h, m, s, ms, us;
}Tempo;

#if defined(_WIN32) || defined(_WIN64)
	#include <windows.h>
    #define SISTEMA "Windows"
    Tempo * tempo(){
		Tempo *tmp = new Tempo;

        SYSTEMTIME st;
	    GetSystemTime(&st);
		tmp->us = 0;
		tmp->ms = st.wMilliseconds;
		tmp->s = st.wSecond;
		tmp->m = st.wMinute;
		tmp->h = st.wHour;

		return tmp;
	}
#elif __linux__
    #include <sys/time.h>
    #define SISTEMA "Linux"
    Tempo * tempo(){
		Tempo *tmp = new Tempo;

		struct timeval st;
        gettimeofday(&st, NULL);
        tmp->us = st.tv_usec;
        tmp->ms = tmp->us/1000;
        tmp->us %= 1000;
        tmp->s = st.tv_sec;
        if (tmp->us < 0){
            tmp->us += 1000;
            tmp->ms += 999;
            tmp->s--;
        }
        tmp->h = tmp->s / 3600;
        tmp->m = (tmp->s / 60) % 60;
        tmp->s %= 60;

		return tmp;
	}
#else
    #error Sistema Operacional nao identificado!!!
#endif

int main(int n_arg, char **args){
    if (n_arg < 2)
    {
        printf("Erro no numero de argumentos, favor informar o nome do arquivo de configuracao de entrada.\n");
        help(args[0]);
        return -1;
    }
#ifdef _OPENMP
    unsigned int NUM_THREADS = 1;
    omp_set_num_threads(NUM_THREADS);
#endif

    unsigned int CONTROLE = 0; // Variavel de controle para a simulacao:
//    bit 0 -> Simulacao nova,
//    bit 1 -> Parada por passo de tempo,
//    bit 2 -> Parada por tempo,
//    bit 3 -> Parada por velocidade,
//    bit 4 -> Parada por energia,
//    bit 5 -> Vibracao da parede.
//    bit 6 -> Nao intruso na simulacao
//    bit 7 -> Executa a SobeIntruso
//    bit 8 -> Empurra a Caixa
//    bit 9 -> Condicao periodica de contorno em x
    unsigned int n_graos = 0;
    double altura = 6.0, largura = 0.15; // Tamanho da caixa de simulacao
    unsigned long long int paradap = 0; // Condicao de parada
    double paradat = 0.0, paradav = 0.0, paradae = 0.0; // Condicoes de parada
    double vibra[2] = {0.0, 0.0}; // Gamma de vibracao e período;
    double velocidadeparede = 0.0; // Velocidade de subida da caixa.

    for (int i = 1; i < n_arg; i++){
        if ((strcmp(args[i], "-h") == 0) ||(strcmp(args[i], "--help") == 0)){
            help(args[0]);
            return 0;
        }
        if (strcmp(args[i], "-n") == 0){ // Numero de Threads associado ao OpenMP, por padrao serial: 1 thread.
#ifdef _OPENMP
            NUM_THREADS = atoi(args[i+1]);
            omp_set_num_threads(NUM_THREADS);
#endif
        }
        if (strcmp(args[i], "-N") == 0){ // Simulacao com nova base de n_graos.
            CONTROLE |= 0x01;
            n_graos = atoi(args[i+1]);
            semente = atoi(args[i+2]);
            altura = atof(args[i+3]);
            largura = atof(args[i+4]);
            largura = ((largura == 0.0) || (altura == 0.0)) ? 0.15 : largura;
            altura = (altura == 0.0) ? 6.0 : altura;
        }
        if (strcmp(args[i], "-p") == 0){ // Parada pelo passo de tempo
            paradap = atoll(args[i+1]);
            CONTROLE |= (paradap != 0) ? 0x02 : 0x00;
        }
        if (strcmp(args[i], "-t") == 0){ // Parada pelo tempo de simulacao
            paradat = atof(args[i+1]);
            CONTROLE |= (paradat != 0.0) ? 0x04 : 0x00;
        }
        if (strcmp(args[i], "-v") == 0){ // Parada pela velocidade do intruso
            paradav = atof(args[i+1]);
            CONTROLE |= (paradav != 0.0) ? 0x08 : 0x00;
        }
        if (strcmp(args[i], "-e") == 0){ // Parada pela energia total
            paradae = atof(args[i+1]);
            CONTROLE |= (paradae != 0.0) ? 0x10 : 0x00;
        }
        if (strcmp(args[i], "-s") == 0){ // Vibracao do sistema
            vibra[0] = atof(args[i+1]); // Aceleração
            vibra[1] = atof(args[i+2]); // Periodo
            CONTROLE |= ((vibra[0] != 0.0) && (vibra[1] != 0.0)) ? 0x20 : 0x00;
        }
        if (strcmp(args[i], "-I") == 0){ // Nao intruso na simulacao
            CONTROLE |= 0x40;
        }
        if (strcmp(args[i], "-T") == 0){ // Sobe o Intruso
            CONTROLE |= 0x80;
        }
        if (strcmp(args[i], "-E") == 0){ // Sobe o Intruso
            velocidadeparede = atof(args[i+1]); // Velocidade de subida
            CONTROLE |= (velocidadeparede != 0.0) ? 0x100 : 0x00;
        }
        if (strcmp(args[i], "-P") == 0){ // Condicao Periodica de contorno em x
            CONTROLE |= 0x200;
        }
    }

    PREDITOR_CORRETOR *Caixa = NULL;
    char tmp[100];

    if (CONTROLE & 0x01){ // Criacao de uma nova caixa e simulacao de 1 segundo de assentamento
        sprintf(tmp, "%s%s", "c", args[1]);
        FILE *Arquivo = fopen(tmp,"w");
        fprintf(Arquivo, "#Sistema: %s\n", SISTEMA);
        fprintf(Arquivo, "#");
        for (int i = 0; i < n_arg; i++)
            fprintf(Arquivo, "%s ", args[i]);
        fprintf(Arquivo, "\n");
        fprintf(Arquivo, "#semente: %d\n", semente);
#ifdef _OPENMP
        fprintf(Arquivo, "#processadores: %d, max_theads: %d, threads: %d\n", omp_get_num_procs(), omp_get_max_threads(), NUM_THREADS);
#endif

        Caixa = new PREDITOR_CORRETOR(largura, altura, n_graos);
        sprintf(tmp, "%s_%s", "inicial", args[1]);
        Caixa->Salvar_Configuracao(tmp);
        double mediae[4] = {0.0, 0.0, 0.0, 0.0}; // Media das energias do sistema no instante de PERIODO passos de tempo

        fprintf(Arquivo, "#dt: %e\n", Caixa->dt);
        printf("dt: %e\n", Caixa->dt);

        clock_t t_process = clock();
        Tempo *T1 = tempo();
        fprintf(Arquivo, "#Inicio: %02d:%02d:%02d:%03d\n", T1->h, T1->m, T1->s, T1->ms);
        printf("Inicio da base: %02d:%02d:%02d:%03d\n", T1->h, T1->m, T1->s, T1->ms);
        delete T1;

        fprintf(Arquivo, "#Tempo Energia(C, R, E, G)\n");

        if (CONTROLE & 0x80) // Condicao de colocar o intruso no topo da pilha granular
            Caixa->SobeIntruso();

        while (Caixa->passos < 10 * PERIODO){ // 10⁴ passos de DM, com dt=1/10 no cotato, indica aproximadamente 10^3 contatos até estabilizar
            if (Caixa->passos % 100 == 0){ // A cada 10 contatos, procura de vizinhos
                Caixa->Procura_Vizinhos2((CONTROLE & 0x200) ? 1 : 0);
            }
            Caixa->Preditor(((CONTROLE & 0x40) ? 0 : 1)+((CONTROLE & 0x200) ? 2 : 0));
            Caixa->Detectar_Contatos((CONTROLE & 0x200) ? 1 : 0);
            Caixa->Calculo_Forca((CONTROLE & 0x200) ? 1 : 0);
            Caixa->Corretor((CONTROLE & 0x40) ? 0 : 1);
            Caixa->passos++;

            mediae[0] += Caixa->Energia_Cinetica();
            mediae[1] += Caixa->Energia_Rotacional();
            mediae[2] += Caixa->Energia_Elastica();
            mediae[3] += Caixa->Energia_Gravitacional();

            if (Caixa->passos % (PERIODO/10) == 0){ // De 10³ em 10³ passos de DM, salvar as configurações
                for (int i = 0; i < 4; i++)
                    mediae[i] /= PERIODO;
                printf("\r%e passos de termalizacao... Tempo de sistema: %e... Altura do intruso: %e\0", Caixa->passos/1.0, Caixa->passos * Caixa->dt, Caixa->Intruso.Linear.Posicao.y);
                fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e\n", Caixa->passos * Caixa->dt, mediae[0], mediae[1], mediae[2], mediae[3]);
                //fflush(stdout);
                //fflush(Arquivo);

                for (int i = 0; i < 4; i++)
                    mediae[i] = 0;

                if (Caixa->passos % (PERIODO * 10) == 0){
                    sprintf(tmp, "%s_termal_%llu.dat", args[1], Caixa->passos);
                    Caixa->Salvar_Configuracao(tmp);
                    T1 = tempo();
                    fprintf(Arquivo, "#Tempo decorrido: %02d:%02d:%02d:%03d\n#Tempo de maquina[s]: %.10e\n", T1->h, T1->m, T1->s, T1->ms, float(clock() -t_process)/CLOCKS_PER_SEC);
                    delete T1;
                }
            }
        }
        Caixa->passos = 0;
        Caixa->Salvar_Configuracao(args[1]);
        fprintf(Arquivo, "#Passos: %llu\n", Caixa->passos);

        T1 = tempo();
        fprintf(Arquivo, "#Fim\n#Tempo real: %02d:%02d:%02d:%03d\n#Tempo de maquina[s]: %.10e\n", T1->h, T1->m, T1->s, T1->ms, float(clock() -t_process)/CLOCKS_PER_SEC);
        printf("\nFim da base: %02d:%02d:%02d:%03d\n", T1->h, T1->m, T1->s, T1->ms);
        delete T1;

        fclose(Arquivo);
    }
    else { // Se for a continuacao de uma simulacao existente, carregar seus parmetros
        FILE *Config;
        unsigned int p = 1, n, t, e = 0, c;
        unsigned long long int params[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        double **dados = NULL, *f_contatos = NULL;
        int **contatos = NULL;

        Config = fopen(args[1],"r");
        {
            char tmpc = 0;
            do {
                tmpc = fgetc(Config);
                e += (tmpc == 32) ? 1 : 0;
            } while(tmpc != 10);
            rewind(Config);
        }

        if (e == 1)
            fscanf(Config, "%d %d\n", &n, &t);
        else if (e >= 2){
            fscanf(Config, "%d", &n);
            fscanf(Config, "%d", &t);
            fscanf(Config, "%d", &p);
            for (int i = 0; i < e -2; i++)
                fscanf(Config, "%llu", &params[i]);
        }

        dados = new double*[t];
        if (p == 1) {
            for (int i = 0; i < t; i++){
                dados[i] = new double[3];
                fscanf(Config, "%le %le %le\n", &dados[i][0], &dados[i][1], &dados[i][2]);
            }
        }
        else {
            // p == 2, indicando apenas as posiçoes, velocidades, acelerações e raios salvos
            for (int i = 0; i < t; i++){
                dados[i] = new double[10];
                fscanf(Config, "%le %le %le %le %le %le %le %le %le %le\n", &dados[i][0], &dados[i][1], &dados[i][2], &dados[i][3], &dados[i][4], &dados[i][5], &dados[i][6], &dados[i][7], &dados[i][8], &dados[i][9]);
            }
            // p == 3, indicando os valores de reação tangencial anteriores
            if (p == 3){
                f_contatos = new double[params[1]];
                contatos = new int*[params[1]];
                for (int i = 0; i < params[1]; i++){
                    contatos[i] = new int[2];
                }
                for (int i = 0; i < params[1]; i++){
                    fscanf(Config, "%i %i %le\n", &contatos[i][0], &contatos[i][1], &f_contatos[i]);
                }
            }
        }
        fclose(Config);
        if (p <= 2)
            Caixa = new PREDITOR_CORRETOR(p, n, dados);
        else
            Caixa = new PREDITOR_CORRETOR(p, n, dados, params[1], contatos, f_contatos);
        if (dados != NULL){
            for (int i = 0; i < t; i++)
                delete [] dados[i];
            delete [] dados;
        }
        if (p == 3){
            if (contatos != NULL){
                for (int i = 0; i < params[1]; i++)
                    delete [] contatos[i];
                delete [] contatos;
            }
            if (f_contatos != NULL){
                delete [] f_contatos;
            }
        }
        Caixa->passos = params[0];
    }

    FILE *Arquivo, *Arquivo2, *Arquivo3;
    char tipo[10];
    sprintf(tipo, (CONTROLE & 0x01) ? "w" : "a+");
    sprintf(tmp, "%s%s", "e", args[1]);
    Arquivo = fopen(tmp, tipo);
    sprintf(tmp, "%s%s", "i", args[1]);
    Arquivo2 = fopen(tmp, tipo);
    sprintf(tmp, "%s%s", "p", args[1]);
    Arquivo3 = fopen(tmp, tipo);
    fprintf(Arquivo, "#Sistema %s\n", SISTEMA);
    fprintf(Arquivo, "#");
    for (int i = 0; i < n_arg; i++)
        fprintf(Arquivo, "%s ", args[i]);
    fprintf(Arquivo, "\n");
    if (CONTROLE & 0x01)
        fprintf(Arquivo, "#semente: %d\n", semente);
#ifdef _OPENMP
    fprintf(Arquivo, "#processadores: %d, max_theads: %d, threads: %d\n", omp_get_num_procs(), omp_get_max_threads(), NUM_THREADS);
#endif

    clock_t t_process = clock();
    Tempo *T1 = tempo();
    fprintf(Arquivo, "#Inicio: %02d:%02d:%02d:%03d\n", T1->h, T1->m, T1->s, T1->ms);
    printf("Inicio: %02d:%02d:%02d:%03d\n", T1->h, T1->m, T1->s, T1->ms);
    delete T1;

    try{
        if (CONTROLE & 0x20){
            // Ajustando o passo de tempo para simulações com dinâmica nos grãos
            Caixa->dt /= 10;
            vibra[1] *= 10;
            fprintf(Arquivo, "#Vibracao da caixa: Aceleração: %e Periodo: %e", vibra[0], vibra[1]);
            // Ajusta a amplitude de vibração do sistema
            vibra[0] *= Caixa->gravidade*vibra[1]*vibra[1]*Caixa->dt*Caixa->dt/(4*PI*PI); // Amplitude
            fprintf(Arquivo, " Amplitude: %e\n", vibra[0]);
            // Posicionando as bordas da caixa
            for (int i = 0; i < 4; i++)
                Caixa->p0y[i] -= vibra[0] * sin(2*PI*(Caixa->passos)/(vibra[1]));
        }
        else
            fprintf(Arquivo, "#Caixa estatica (sem vibracoes)\n");

        fprintf(Arquivo, "#dt: %e\n", Caixa->dt);
        fprintf(Arquivo, "#Massa: %e Raio: %e\n", Caixa->Intruso.massa, Caixa->Intruso.Geometria.raio);
        printf("dt: %e numero_graos:%d\n", Caixa->dt, Caixa->numero_graos);
        fprintf(Arquivo, "#Tempo Energia(C, R, E, G)\n");
        fprintf(Arquivo2, "#Tempo Posicao(x,y,theta) Velocidade(x,y,omega) Aceleracao(x,y,alpha)\n");
        fprintf(Arquivo3, "#Tempo Altura Interpenetracao Interpenetracao_Total Compactacao Coordenacao Coordenacao_Total\n");

        double mediae[4] = {0.0, 0.0, 0.0, 0.0}; // Media das energias do sistema no instante de PERIODO passos de tempo

        unsigned long long int PARADAe = 0, PARADAv = 0; // Contadores da condicao de parada: energia total e velocidade do intruso.

        if (CONTROLE & 0x100){ // Subir parede
            double maxx = Caixa->Intruso.Linear.Posicao.x +Caixa->Intruso.Geometria.raio, minx = Caixa->Intruso.Linear.Posicao.x -Caixa->Intruso.Geometria.raio;
            double maxy = Caixa->Intruso.Linear.Posicao.y +Caixa->Intruso.Geometria.raio, miny = Caixa->Intruso.Linear.Posicao.y -Caixa->Intruso.Geometria.raio;
            for(int i = 0; i < Caixa->numero_graos; i++){
                maxx = maxx < Caixa->Graos[i].Linear.Posicao.x +Caixa->Graos[i].Geometria.raio ? Caixa->Graos[i].Linear.Posicao.x +Caixa->Graos[i].Geometria.raio : maxx;
                minx = minx > Caixa->Graos[i].Linear.Posicao.x -Caixa->Graos[i].Geometria.raio ? Caixa->Graos[i].Linear.Posicao.x -Caixa->Graos[i].Geometria.raio : minx;
                maxy = maxy < Caixa->Graos[i].Linear.Posicao.y +Caixa->Graos[i].Geometria.raio ? Caixa->Graos[i].Linear.Posicao.y +Caixa->Graos[i].Geometria.raio : maxy;
                miny = miny > Caixa->Graos[i].Linear.Posicao.y -Caixa->Graos[i].Geometria.raio ? Caixa->Graos[i].Linear.Posicao.y -Caixa->Graos[i].Geometria.raio : miny;
            }
            minx -= Caixa->Bordas[0].Geometria.raio*1;
            miny -= Caixa->Bordas[0].Geometria.raio*1;
            maxx += Caixa->Bordas[0].Geometria.raio*1;
            maxy += Caixa->Bordas[0].Geometria.raio*1;

            // Definindo as posicoes das paredes
            Caixa->Bordas[0].Linear.Posicao.x = minx;
            Caixa->Bordas[0].Linear.Posicao.y = miny;
            Caixa->Bordas[1].Linear.Posicao.x = maxx;
            Caixa->Bordas[1].Linear.Posicao.y = miny;
            Caixa->Bordas[2].Linear.Posicao.x = minx;
            Caixa->Bordas[2].Linear.Posicao.y = maxy;
            Caixa->Bordas[3].Linear.Posicao.x = maxx;
            Caixa->Bordas[3].Linear.Posicao.y = maxy;

            Caixa->gravidade = 0.0;
        }

        if (CONTROLE & 0x80) // Condicao de colocar o intruso no topo da pilha granular
            Caixa->SobeIntruso();

	    while (((CONTROLE & 0x02)&&(Caixa->passos < paradap))||((CONTROLE & 0x04)&&(Caixa->passos * Caixa->dt < paradat))||((CONTROLE & 0x08)&&(PARADAv < PERIODO))||((CONTROLE & 0x10)&&(PARADAe < PERIODO))){ // Condicoes de parada estabelecidads.
            if (Caixa->passos % 100 == 0)
                Caixa->Procura_Vizinhos2((CONTROLE & 0x200) ? 0 : 1);
            if (CONTROLE & 0x20) // Condicao de vibracao da parede.
                Caixa->Vibra_Parede(vibra[0], vibra[1]);
            if (CONTROLE & 0x100) // Subir parede
                Caixa->Subir_Parede(velocidadeparede);
            Caixa->Preditor(((CONTROLE & 0x40) ? 0 : 1) +(CONTROLE & 0x200) ? 2 : 0); //Passando ou nao o intruso para a simulacao
            Caixa->Detectar_Contatos((CONTROLE & 0x200) ? 1 : 0);
            Caixa->Calculo_Forca((CONTROLE & 0x200) ? 1 : 0);
            Caixa->Corretor((CONTROLE & 0x40) ? 0 : 1); //Passando ou nao o intruso para a simulacao
            Caixa->passos++;

            PARADAv = (fabs(Caixa->Intruso.Linear.Velocidade.norma()) < paradav) ? PARADAv + 1 : 0;
            PARADAe = ((Caixa->Energia_Cinetica()+Caixa->Energia_Rotacional()+Caixa->Energia_Elastica()+Caixa->Energia_Gravitacional()) < paradae) ? PARADAe + 1 : 0;

//FORMAS LEGIVEIS DAS CONDICOES DE PARADA
/*            if ((CONTROLE & 0x02)&&(Caixa->passos > paradap))
                break;
            if ((CONTROLE & 0x04)&&(Caixa->passos * Caixa->dt > paradat))
                break;
            if ((CONTROLE & 0x08)&&(PARADAv > PERIODO))
                break;
            if ((CONTROLE & 0x10)&&(PARADAe > PERIODO))
                break;*/

		    mediae[0] += Caixa->Energia_Cinetica();
            mediae[1] += Caixa->Energia_Rotacional();
            mediae[2] += Caixa->Energia_Elastica();
            mediae[3] += Caixa->Energia_Gravitacional();

            if (Caixa->passos % PERIODO == 0){
                for (int i = 0; i < 4; i++)
                    mediae[i] /= PERIODO;
                printf("\r%e passos... Tempo de sistema: %e...\0", Caixa->passos/1.0, Caixa->passos * Caixa->dt);
                fprintf(Arquivo, "%.10e %.10e %.10e %.10e %.10e\n", Caixa->passos * Caixa->dt, mediae[0], mediae[1], mediae[2], mediae[3]);
                fprintf(Arquivo2, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Caixa->passos * Caixa->dt, Caixa->Intruso.Linear.Posicao.x, Caixa->Intruso.Linear.Posicao.y, Caixa->Intruso.Angular.Posicao.x, Caixa->Intruso.Linear.Velocidade.x, Caixa->Intruso.Linear.Velocidade.y, Caixa->Intruso.Angular.Velocidade.x, Caixa->Intruso.Linear.Aceleracao.x, Caixa->Intruso.Linear.Aceleracao.y, Caixa->Intruso.Angular.Aceleracao.x);
                fprintf(Arquivo3, "%.10e %.10e %.10e %.10e %.10e %.10e %.10e\n", Caixa->passos * Caixa->dt, Caixa->Intruso.Linear.Posicao.y, Caixa->Interpenetracao(), Caixa->Area_Interpenetracao_Total(), Caixa->Compactacao(), Caixa->Coordenacao(0), Caixa->Coordenacao(1));

                for (int i = 0; i < 4; i++)
                    mediae[i] = 0;

                if (Caixa->passos % (PERIODO * 10) == 0){
                    sprintf(tmp, "%s_simul_%010llu.dat", args[1], Caixa->passos);
                    Caixa->Salvar_Configuracao(tmp);
                    T1 = tempo();
                    fprintf(Arquivo, "#Tempo decorrido: %02d:%02d:%02d:%03d\n#Tempo de maquina[s]: %.10e\n", T1->h, T1->m, T1->s, T1->ms, float(clock() -t_process)/CLOCKS_PER_SEC);
                    delete T1;
                }

                //fflush(stdout);
                //fflush(Arquivo);
                //fflush(Arquivo2);
                //fflush(Arquivo3);
			}
        }
    }
    catch(...){
        fprintf(Arquivo, "#Erro - Tempo real: %02d:%02d:%02d:%03d\n#Tempo de maquina[s]: %.10e\n", T1->h, T1->m, T1->s, T1->ms, float(clock() -t_process)/CLOCKS_PER_SEC);
        fprintf(stderr, "ERRO INESPERADO!\n");
    }

    sprintf(tmp, "%s_%s", "final", args[1]);
    Caixa->Salvar_Configuracao(tmp);
    fprintf(Arquivo, "#Passos: %llu\n", Caixa->passos);

	T1 = tempo();
    fprintf(Arquivo, "#Fim\n#Tempo real: %02d:%02d:%02d:%03d\n#Tempo de maquina[s]: %.10e\n", T1->h, T1->m, T1->s, T1->ms, float(clock() -t_process)/CLOCKS_PER_SEC);
    printf("Fim: %02d:%02d:%02d:%03d\n", T1->h, T1->m, T1->s, T1->ms);
	delete T1;

    fclose(Arquivo);
    fclose(Arquivo2);
    fclose(Arquivo3);

    return 0;
}

void help(const char * nome){
    printf("%s CONFIGURACAO [OPCOES]\n\n", nome);
    printf("Simulador de Sistemas Granulares confinados para estudo de um intruso no meio granular.\n");
    printf("O programa deve ser chamado com pelo menos um argumento, que e o nome da configuracao inicial de simulacao.\n");
    printf("Parmetros de passagem adicionais para a simulacao:\n");
    printf("    -h :            Abre esta ajuda.\n");
    printf("    -n []:          Numero de threads disparadas nas clausulas OpenMP (Numero natural).\n");
    printf("    -N [] [] [] []: Gera uma configuracao com o numero de graos (Numero natural) a partir de uma configuracao preparada aleatoriamente com base inicial de numeros aleatorios (Numero natural) de altura (Numero real) e largura (Numero real).\n");
    printf("    -p []:          Criterio de parada por passo de tempo (Numero natural).\n");
    printf("    -t []:          Criterio de parada por equivalente de tempo de simulacao (Numero real).\n");
    printf("    -v []:          Criterio de parada por velocidade (Numero real) \"estacionaria\" do intruso em %i passos de tempo.\n", PERIODO);
    printf("    -e []:          Criterio de parada por energia (Numero real) \"estacionaria\" do sistema em %i passos de tempo.\n", PERIODO);
    printf("    -s [] []:       Vibra a parede com aceleração em função da gravidade (Numero real) e periodo baseado no periodo natural de contato granular (Numero real).\n");
    printf("    -E []:          Empurra a parede com uma velocidade imposta (Numero real).\n");
    printf("    -I :            Fixa a posicao do intruso.\n");
    printf("    -T :            Coloca a posicao do intruso no topo da rede granular.\n");
    printf("    -P :            Coloca condicao periodica de contorno no eixo x.\n");
    printf("\n");
}
