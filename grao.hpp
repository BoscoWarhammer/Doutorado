#define N_CONTATOS_MAX 100
#define N_VIZINHOS_MAX 1000

class GEOMETRIA;
class LINEAR;
class ANGULAR;
class VETOR;
class GRAO;
struct nt;

class GEOMETRIA {
    public:
        double raio;               // Circular em forma de disco
    public:
        GEOMETRIA ();             // Construtor padrao do objeto
        ~GEOMETRIA();             // Destrutor padrao do objeto
        GEOMETRIA (const GEOMETRIA&);    // Construtor do objeto de definicao de geometria do grao em relacao `a outro GRAO
        GEOMETRIA (double);        // Construtor do objeto de definicao de geometria do grao em relacao de disco bi-dimensional ou esfera
};

class VETOR {
    public:
        double x, y;              // Atributos de um vetor tridimencional
    public:
        VETOR ();                 // Construtor padrao do objeto
        ~VETOR ();                // Destrutor padrao do objeto
        VETOR (const VETOR&);            // Construtor do objeto em relacao `a outro VETOR
        VETOR (double, double);    // Construtor do objeto em relacao `a 3 numeros
        VETOR (double *);           // Construtor do objeto em relacao `a outro vetor
        VETOR operator + (VETOR); // Definindo operacao de soma de VETOR em relacao ao outro VETOR
        VETOR operator - (VETOR); // Definindo operacao de subtracao de VETOR em relacao ao outro VETOR
        VETOR operator = (const VETOR&); // Definindo operacao de atribuicao de VETOR em relacao ao outro VETOR
        VETOR operator + (double*); // Definindo operacao de soma de VETOR em relacao `a um vetor
        VETOR operator - (double*); // Definindo operacao de subracao de VETOR em relacao `a um vetor
        VETOR operator = (double*); // Definindo operacao de atribuicao de VETOR em relacao `a um vetor
        VETOR operator += (VETOR); // Definindo operacao de atribuicao somada de um VETOR em relacao `a um vetor
        VETOR operator -= (VETOR); // Definindo operacao de atribuicao subtraida de um VETOR em relacao `a um vetor
        VETOR operator += (double *); // Definindo operacao de atribuicao somada de um VETOR em relacao `a um vetor
        VETOR operator -= (double *); // Definindo operacao de atribuicao subtraida de um VETOR em relacao `a um vetor
        double norma();
};

class VETOR_POLAR : public VETOR {
    public:
        double r, t;
    public:
        VETOR_POLAR(double, double);
};

class LINEAR {
    public:
        VETOR Posicao;            // Atributos espaciais
        VETOR Velocidade;         // Atributos derivadas primeiras temporais
        VETOR Aceleracao;         // Atributos derivadas segundas temporais
    public:
        LINEAR();                 // Construtor padrao do objeto
        LINEAR(VETOR, VETOR, VETOR); // Construtor do objeto definindo os vetores posicao, velocidade e aceleracao, respectivamente
        ~LINEAR();                // Destrutor padrao do objeto
};

class ANGULAR{
    public:
        VETOR Posicao;            // Atributos espaciais de vetor normal ao plano de rotacao
        VETOR Velocidade;         // Atributos derivadas primeiras temporais de vetor normal ao plano de rotacao
        VETOR Aceleracao;         // Atributos derivadas segundas temporais de vetor normal ao plano de rotacao
    public:
        ANGULAR();                // Construtor padrao do objeto
        ANGULAR(VETOR, VETOR, VETOR); // Construtor do objeto definindo os vetores posicao, velocidade e aceleracao, respectivamente
        ~ANGULAR();               // Destrutor padrao do objeto
};

struct nt{
    public:
        double normal;             // Coeficiente Normal
        double tangencial;         // Coeficiente Tangencial
};

struct REACAO_TANGENCIAL{
    public:
        double *forca;   // Modulo da forca de contato no passo de tempo anterior
        int *contato;      // Lista do grao de contato no passo de tempo anterior
};

class GRAO {
    public:
        GEOMETRIA Geometria;      // Objeto que descreve a geometria do Grao
        LINEAR Linear, PrevistoL;            // Objeto que descreve a posicao, velocidade e aceleracao lineares do grao e suas previsoes
        ANGULAR Angular, PrevistoA;          // Objeto que descreve a posicao, velocidade e aceleracao angulares do grao e suas previsoes
		REACAO_TANGENCIAL Reacao_Tangencial; //Estrutura que descreve o historico da reacao tangencial da particula no passo de tempo passado
        VETOR Forca; //Previsao das forcas utilizadas no corretor
		double Torque; //Previsao dos torques utilizados no corretor -- Um numero real pois o sistema e em duas dimensoes, portanto o torque so pode estar em uma unica dimensao
        double *atrito;            // Atributo que relaciona as constantes de atrito entre as superficies dos graos com os meios de contato por meio de matriz
        nt restituicao;           // Estrutura que tem a funcao de constante de restituicao relacionada ao grao
        double massa;              // Atributo que contem a massa do grao
        double densidade;          // Atributo que contem a densidade do grao
        double inercia;            // Atributo que contem a inercia do grao
        VETOR CentroMassa;        // Objeto que contem a regiao geometrica do Centro de Massa do Grao
        VETOR CentroGeometrico;   // Objeto que contem a regiao de centro de geometria do grao
        nt amortecimento;         // Estrutura que contem a constante de amortecimento da equacao diferencial
        nt mola;                   // Estrutura que contem a constante de mola da equacao diferencial
        int *vizinhos;             // Atributo que contem no indice 0 o numero de vizinhos e em seguida a identificacao dos vizinhos
        int *contatos;             // Atributo que contem no indice 0 o numero de contatos e em seguida a identificacao dos contatos
    public:
        GRAO ();                  // Construtor padrao do objeto
        ~GRAO ();                 // Destrutor padrao do objeto
        double Distancia (const GRAO&);   // Distancia entre dois graos
        double Interpenetracao (const GRAO&); // Vetor de interpenetracao entre os graos, este e um outro de referencia
};
