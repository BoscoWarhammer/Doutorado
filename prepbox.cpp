#include "grao.cpp"

int semente = 0;

double random(int *base);

class CAIXA {
    public:
        GRAO Bordas[4]; // Bordas definidas inicialmente apartir do intruso
        bool caixa_pronta;
        double largura, altura;
        GRAO Intruso;
        GRAO *Graos;
        unsigned int numero_graos;
        double densidade;
    public:
        CAIXA ();     // Construtor padrao
        CAIXA (double, double, unsigned int);     // Construtor definindo a largura e altura
        CAIXA (int, int, double**);
        CAIXA (int, int, double**, int, int**, double*);
        ~CAIXA ();
        double Densidade ();
        int GerarCaixa ();
        int SobeIntruso();
        int InserirIntruso (GRAO);
        int InserirIntruso (double, VETOR); // Insere o grao apartir de um raio e uma posicao
        int InserirGraos (GRAO *);
        double Condicao_Periodica_x(double);
        GRAO* Preencher_Caixa();
};

CAIXA::CAIXA(int p, int n, double **dados, int numero_contatos, int **contatos, double *f_contatos) : CAIXA(p, n, dados){
    Graos = NULL;
    numero_graos = n;
    Graos = new GRAO[numero_graos];
    if (Graos == NULL){
        printf("Erro de alocacao de memoria, sem espaco para alocar os graos.\n");
        exit(EXIT_FAILURE);
    }

     if (p == 3){
        for (int i = 0; i < numero_graos; i++){
            Graos[i].contatos[0] = 0;
        }

        for (int i = 0; i < numero_contatos; i++){
            if (contatos[i][0] != -1){
                Graos[contatos[i][0]].contatos[++Graos[contatos[i][0]].contatos[0]] = contatos[i][1];
                Graos[contatos[i][0]].Reacao_Tangencial.forca[Graos[contatos[i][0]].contatos[0]] = f_contatos[i];
            } else {
                Intruso.Reacao_Tangencial.contato[++Intruso.contatos[0]] = contatos[i][1];
                Intruso.Reacao_Tangencial.forca[Intruso.contatos[0]] = f_contatos[i];
            }
        }
    }

    Intruso.Linear.Posicao.x = dados[numero_graos][0];
    Intruso.Linear.Posicao.y = dados[numero_graos][1];
    if (p == 1)
        Intruso.Geometria.raio = dados[numero_graos][2];
    else if (p >= 2){
        Intruso.Angular.Posicao.x = dados[numero_graos][2];
        Intruso.Linear.Velocidade.x = dados[numero_graos][3];
        Intruso.Linear.Velocidade.y = dados[numero_graos][4];
        Intruso.Angular.Velocidade.x = dados[numero_graos][5];
        Intruso.Linear.Aceleracao.x = dados[numero_graos][6];
        Intruso.Linear.Aceleracao.y = dados[numero_graos][7];
        Intruso.Angular.Aceleracao.x = dados[numero_graos][8];
        Intruso.Geometria.raio = dados[numero_graos][9];
    }
    Intruso.massa = Intruso.Geometria.raio * Intruso.Geometria.raio;
    Intruso.inercia = (Intruso.massa * Intruso.Geometria.raio * Intruso.Geometria.raio) / 2;
    Intruso.mola.normal = 1000;
    Intruso.mola.tangencial = 0.75 * Intruso.mola.normal;
    Intruso.amortecimento.normal = 2 * sqrt(Intruso.mola.normal * Intruso.massa);
    Intruso.amortecimento.tangencial = 0.75 * Intruso.amortecimento.normal;

    double maior_altura = dados[numero_graos +1][1], RAIO = 0;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = 0; i < numero_graos; i++){
        Graos[i].Linear.Posicao.x = dados[i][0];
        Graos[i].Linear.Posicao.y = dados[i][1];
        if (p == 1)
            Graos[i].Geometria.raio = dados[i][2];
        else if (p >= 2){
            Graos[i].Angular.Posicao.x = dados[i][2];
            Graos[i].Linear.Velocidade.x = dados[i][3];
            Graos[i].Linear.Velocidade.y = dados[i][4];
            Graos[i].Angular.Velocidade.x = dados[i][5];
            Graos[i].Linear.Aceleracao.x = dados[i][6];
            Graos[i].Linear.Aceleracao.y = dados[i][7];
            Graos[i].Angular.Aceleracao.x = dados[i][8];
            Graos[i].Geometria.raio = dados[i][9];
        }
        Graos[i].massa = Graos[i].Geometria.raio * Graos[i].Geometria.raio;
        Graos[i].inercia = (Graos[i].massa * Graos[i].Geometria.raio * Graos[i].Geometria.raio) / 2;
        Graos[i].mola.normal = 1000;
        Graos[i].mola.tangencial = 0.75 * Graos[i].mola.normal;
        Graos[i].amortecimento.normal = 2 * sqrt(Graos[i].mola.normal * Graos[i].massa);
        Graos[i].amortecimento.tangencial = 0.75 * Graos[i].amortecimento.normal;
    }

    for(int i = 0; i < 4; i++){
        Bordas[i].Linear.Posicao.x = dados[numero_graos +i+1][0];
        Bordas[i].Linear.Posicao.y = dados[numero_graos +i+1][1];
        if (p == 1)
            Bordas[i].Geometria.raio = dados[numero_graos +i+1][2];
        else if (p >= 2)
            Bordas[i].Geometria.raio = dados[numero_graos +i+1][9];
        Bordas[i].massa = Bordas[i].Geometria.raio * Bordas[i].Geometria.raio;
        Bordas[i].inercia = Bordas[i].massa * Bordas[i].Geometria.raio * Bordas[i].Geometria.raio / 2;
        Bordas[i].mola.normal = 1000;
        Bordas[i].mola.tangencial = 0.75 * Bordas[i].mola.normal;
        Bordas[i].amortecimento.normal = 2 * sqrt(Bordas[i].mola.normal * Bordas[i].massa);
    }
}

CAIXA::CAIXA(int p, int n, double **dados){
    Graos = NULL;
    numero_graos = n;
    Graos = new GRAO[numero_graos];
    if (Graos == NULL){
        printf("Erro de alocacao de memoria, sem espaco para alocar os graos.\n");
        exit(EXIT_FAILURE);
    }

    Intruso.Linear.Posicao.x = dados[numero_graos][0];
    Intruso.Linear.Posicao.y = dados[numero_graos][1];
    if (p == 1)
        Intruso.Geometria.raio = dados[numero_graos][2];
    else if (p >= 2){
        Intruso.Angular.Posicao.x = dados[numero_graos][2];
        Intruso.Linear.Velocidade.x = dados[numero_graos][3];
        Intruso.Linear.Velocidade.y = dados[numero_graos][4];
        Intruso.Angular.Velocidade.x = dados[numero_graos][5];
        Intruso.Linear.Aceleracao.x = dados[numero_graos][6];
        Intruso.Linear.Aceleracao.y = dados[numero_graos][7];
        Intruso.Angular.Aceleracao.x = dados[numero_graos][8];
        Intruso.Geometria.raio = dados[numero_graos][9];
    }
    Intruso.massa = Intruso.Geometria.raio * Intruso.Geometria.raio;
    Intruso.inercia = (Intruso.massa * Intruso.Geometria.raio * Intruso.Geometria.raio) / 2;
    Intruso.mola.normal = 1000;
    Intruso.mola.tangencial = 0.75 * Intruso.mola.normal;
    Intruso.amortecimento.normal = 2 * sqrt(Intruso.mola.normal * Intruso.massa);
    Intruso.amortecimento.tangencial = 0.75 * Intruso.amortecimento.normal;

    double maior_altura = dados[numero_graos +1][1], RAIO = 0;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = 0; i < numero_graos; i++){
        Graos[i].Linear.Posicao.x = dados[i][0];
        Graos[i].Linear.Posicao.y = dados[i][1];
        if (p == 1)
            Graos[i].Geometria.raio = dados[i][2];
        else if (p >= 2){
            Graos[i].Angular.Posicao.x = dados[i][2];
            Graos[i].Linear.Velocidade.x = dados[i][3];
            Graos[i].Linear.Velocidade.y = dados[i][4];
            Graos[i].Angular.Velocidade.x = dados[i][5];
            Graos[i].Linear.Aceleracao.x = dados[i][6];
            Graos[i].Linear.Aceleracao.y = dados[i][7];
            Graos[i].Angular.Aceleracao.x = dados[i][8];
            Graos[i].Geometria.raio = dados[i][9];
        }
        Graos[i].massa = Graos[i].Geometria.raio * Graos[i].Geometria.raio;
        Graos[i].inercia = (Graos[i].massa * Graos[i].Geometria.raio * Graos[i].Geometria.raio) / 2;
        Graos[i].mola.normal = 1000;
        Graos[i].mola.tangencial = 0.75 * Graos[i].mola.normal;
        Graos[i].amortecimento.normal = 2 * sqrt(Graos[i].mola.normal * Graos[i].massa);
        Graos[i].amortecimento.tangencial = 0.75 * Graos[i].amortecimento.normal;
    }

    for(int i = 0; i < 4; i++){
        Bordas[i].Linear.Posicao.x = dados[numero_graos +i+1][0];
        Bordas[i].Linear.Posicao.y = dados[numero_graos +i+1][1];
        if (p == 1)
            Bordas[i].Geometria.raio = dados[numero_graos +i+1][2];
        else if (p >= 2)
            Bordas[i].Geometria.raio = dados[numero_graos +i+1][9];
        Bordas[i].massa = Bordas[i].Geometria.raio * Bordas[i].Geometria.raio;
        Bordas[i].inercia = Bordas[i].massa * Bordas[i].Geometria.raio * Bordas[i].Geometria.raio / 2;
        Bordas[i].mola.normal = 1000;
        Bordas[i].mola.tangencial = 0.75 * Bordas[i].mola.normal;
        Bordas[i].amortecimento.normal = 2 * sqrt(Bordas[i].mola.normal * Bordas[i].massa);
    }
}

CAIXA::CAIXA(){
    Graos = NULL;
    caixa_pronta = false;
    numero_graos = 1000;

    GerarCaixa();
}

CAIXA::CAIXA(double l, double a, unsigned int n){
    Graos = NULL;
    caixa_pronta = false;
    numero_graos = n;
    largura = l;
    altura = a;

    GerarCaixa();
}

double CAIXA::Condicao_Periodica_x(double posicao){
    if (posicao > largura/2){
        return posicao -largura;
    } else if (posicao < -largura/2){
        return posicao +largura;
    } else return posicao;
}

int CAIXA::SobeIntruso(){
    double maior_altura = Graos[0].Linear.Posicao.y, RAIO = 0;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for(int i = 1; i < numero_graos; i++){
        maior_altura = (maior_altura < Graos[i].Linear.Posicao.y) ? Graos[i].Linear.Posicao.y : maior_altura;
        RAIO = (RAIO < Graos[i].Geometria.raio) ? Graos[i].Geometria.raio: RAIO;
    }

    Intruso.Linear.Posicao.y = Intruso.Geometria.raio +maior_altura + 1.5*RAIO;

    return 0;
}

// Rotina que prepara a caixa retangular com os parametros de altura, largura, Intruso, Graos definidos. Erro de colisao no processo resulta em retorno 1.
int CAIXA::GerarCaixa (){
    if ((numero_graos == 0) || (Graos == NULL))
        if (Preencher_Caixa() == NULL)
        {
            printf("Erro! Caixa nao preenchida!\n");
            return -1; // Retorno de caixa sem graos - Provavelmente problema de memoria
        }

    VETOR Min(Intruso.Linear.Posicao.x - Intruso.Geometria.raio, Intruso.Linear.Posicao.y - Intruso.Geometria.raio), Max(Intruso.Linear.Posicao.x + Intruso.Geometria.raio, Intruso.Linear.Posicao.y + Intruso.Geometria.raio + altura);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numero_graos; i++){
        Max.x = (Max.x < Graos[i].Linear.Posicao.x) ? Graos[i].Linear.Posicao.x : Max.x;
        Max.y = (Max.y < Graos[i].Linear.Posicao.y) ? Graos[i].Linear.Posicao.y : Max.y;
        Min.x = (Min.x > Graos[i].Linear.Posicao.x) ? Graos[i].Linear.Posicao.x : Min.x;
        Min.y = (Min.y > Graos[i].Linear.Posicao.y) ? Graos[i].Linear.Posicao.y : Min.y;
    }

    for (int i = 0; i < 4; i++){
        Bordas[i].Geometria.raio = Intruso.Geometria.raio;
        Bordas[i].massa = Bordas[i].Geometria.raio * Bordas[i].Geometria.raio;
        Bordas[i].inercia = Bordas[i].massa * Bordas[i].Geometria.raio * Bordas[i].Geometria.raio / 2;
        Bordas[i].mola.normal = 1000;
        Bordas[i].mola.tangencial = 0.75 * Bordas[i].mola.normal;
        Bordas[i].amortecimento.normal = 2 * sqrt(Bordas[i].mola.normal * Bordas[i].massa);
    }

    Min.x -= Bordas[0].Geometria.raio + (Intruso.Geometria.raio / 3);
    Min.y -= Bordas[0].Geometria.raio + (Intruso.Geometria.raio / 3);
    Max.x += Bordas[0].Geometria.raio + (Intruso.Geometria.raio / 3);
    Max.y += Bordas[0].Geometria.raio + (Intruso.Geometria.raio / 3);

// Definindo as posicoes das paredes
    Bordas[0].Linear.Posicao = Min;
    Bordas[1].Linear.Posicao = Min;
    Bordas[1].Linear.Posicao.x = Max.x;
    Bordas[2].Linear.Posicao = Min;
    Bordas[2].Linear.Posicao.y = Max.y;
    Bordas[3].Linear.Posicao = Max;

    double d = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+: d)
#endif
    for(int i = 0; i < numero_graos; i++)
        d += PI*pow(Graos[i].Geometria.raio,2);
    d += PI*pow(Intruso.Geometria.raio,2);
    densidade = d/((Bordas[0].Linear.Posicao.x - Bordas[1].Linear.Posicao.x - Bordas[0].Geometria.raio - Bordas[1].Geometria.raio) * (Bordas[0].Linear.Posicao.y - Bordas[2].Linear.Posicao.y - Bordas[0].Geometria.raio - Bordas[2].Geometria.raio));

    caixa_pronta = true;
    return 0;
}

GRAO * CAIXA::Preencher_Caixa (){
    double RAIO = 0.0;
    if (Intruso.Geometria.raio == 1.0){
        RAIO = 0.5;
        Intruso.Geometria.raio = 5.0*RAIO;
        Intruso.massa = Intruso.Geometria.raio * Intruso.Geometria.raio;
        Intruso.inercia = (Intruso.massa * Intruso.Geometria.raio * Intruso.Geometria.raio) / 2;
        Intruso.mola.normal = 1000;
        Intruso.mola.tangencial = 0.75 * Intruso.mola.normal;
        Intruso.amortecimento.normal = 2 * sqrt(Intruso.mola.normal * Intruso.massa);
        Intruso.amortecimento.tangencial = 0.75 * Intruso.amortecimento.normal;
        Intruso.Linear.Posicao.x = 0.0;
        Intruso.Linear.Posicao.y = (RAIO +Intruso.Geometria.raio);
    }

    Graos = new GRAO[numero_graos];

    if (Graos == NULL){
        printf("Erro de alocacao de espaco para o numero de graos.\n");
        return NULL;
    }
#ifdef _OPENMP
    #pragma omp parallel for ordered
#endif
    for(int i = 0; i < numero_graos; i++){
		double r = 0.0;
#ifdef _OPENMP
            #pragma omp ordered
#endif
            r = random(&semente);
            Graos[i].Geometria.raio = RAIO * (1 + (r -0.5)*0.025);
            Graos[i].Linear.Posicao.x = 2.5*RAIO * ((i % ((int) (largura / (2*RAIO)))) -((int) (largura / (2*RAIO))/2));
            if (i < numero_graos/2){
                Graos[i].Linear.Posicao.y = Intruso.Linear.Posicao.y + -2.5*RAIO * (i / ((int) (largura / (2*RAIO)))) - RAIO - Intruso.Geometria.raio;
            } else {
                Graos[i].Linear.Posicao.y = 2.5*RAIO * ((i- numero_graos/2) / ((int) (largura / (2*RAIO)))) + RAIO + Intruso.Linear.Posicao.y + Intruso.Geometria.raio;
            }
            Graos[i].massa = Graos[i].Geometria.raio * Graos[i].Geometria.raio;
            Graos[i].inercia = (Graos[i].massa * Graos[i].Geometria.raio * Graos[i].Geometria.raio) / 2;
            Graos[i].mola.normal = 1000;
            Graos[i].mola.tangencial = 0.75 * Graos[i].mola.normal;
            Graos[i].amortecimento.normal = 2 * sqrt(Graos[i].mola.normal * Graos[i].massa);
            Graos[i].amortecimento.tangencial = 0.75 * Graos[i].amortecimento.normal;
    }

    return Graos;
}

CAIXA::~CAIXA(){
    if (Graos != NULL)
        delete[] Graos;
}

double random(int *base){
    int a, b, m;
    a = 843314861;
    b = 453816693;
    m = 1073741824;

    (*base) = (*base) * a + b;

    if (*base < 0) *base = (*base + m) +m;

    return (double (*base) / (2.0* double (m)));
}
