#include "grao.hpp"
#include <math.h>
#include <stdlib.h>

#define NULL __null

const double PI = acos(-1);

GRAO::GRAO(){
    contatos = NULL;
    vizinhos = NULL;
    atrito = NULL;
    Reacao_Tangencial.forca = NULL;
    Reacao_Tangencial.contato = NULL;

    contatos = new int[N_VIZINHOS_MAX];
    vizinhos = new int[N_CONTATOS_MAX];
    Reacao_Tangencial.forca = new double[N_CONTATOS_MAX];
    Reacao_Tangencial.contato = new int[N_CONTATOS_MAX];

    contatos[0] = 0;
    vizinhos[0] = 0;
	Reacao_Tangencial.forca[0] = 0.0;
    Reacao_Tangencial.contato[0] = 0;

    for (int i = 1; i < N_CONTATOS_MAX; i++){
        Reacao_Tangencial.forca[i] = 0.0;
        Reacao_Tangencial.contato[i] = -10;
    }

    Linear.Posicao = {0.0, 0.0};
    Linear.Velocidade = {0.0, 0.0};
    Linear.Aceleracao = {0.0, 0.0};

    Angular.Posicao = {0.0, 0.0};
    Angular.Velocidade = {0.0, 0.0};
    Angular.Aceleracao = {0.0, 0.0};

	Forca = {0.0, 0.0};
	Torque = 0;

    atrito = new double[3];
    atrito[0] = 0.175; // Atrito entre isopor-isopor dividido por 2
    atrito[1] = 0.25; // Atrito entre isopor-aco dividido por 2
    atrito[2] = 0.4; // Atrito entre aco-aco dividido por 2

    densidade = 1/(PI * Geometria.raio * Geometria.raio); // Densidade de 14kg/m^3
    massa = Geometria.raio * Geometria.raio;
    inercia = (massa * Geometria.raio * Geometria.raio) / 2;

    CentroMassa = {0.0, 0.0};
    CentroGeometrico = {0.0, 0.0};

    mola.normal = 1000;
    mola.tangencial = 0.75 * mola.normal;

    amortecimento.normal = 2 * sqrt(mola.normal * massa);
    amortecimento.tangencial = 0.75 * amortecimento.normal;
}

double GRAO::Distancia(const GRAO &a){
    return sqrt(pow(this->Linear.Posicao.x -a.Linear.Posicao.x, 2) +pow(this->Linear.Posicao.y -a.Linear.Posicao.y, 2));
}

double GRAO::Interpenetracao(const GRAO &a){
    if (Distancia(a) < this->Geometria.raio + a.Geometria.raio)
        return (Geometria.raio + a.Geometria.raio -Distancia(a));
    return 0;
}

GEOMETRIA::GEOMETRIA(){
    raio = 1;
}

GEOMETRIA::GEOMETRIA(double r){
    raio = r;
}

GEOMETRIA::GEOMETRIA(const GEOMETRIA &padrao){
    raio = padrao.raio;
}

VETOR::VETOR(){
    x = 0;
    y = 0;
}

VETOR::VETOR(double * componentes){
    this->x = componentes[0];
    this->y = componentes[1];
}

VETOR::VETOR(double x, double y){
    this->x = x;
    this->y = y;
}

VETOR::VETOR(const VETOR &padrao){
    this->x = padrao.x;
    this->y = padrao.y;
}

VETOR VETOR::operator= (const VETOR &padrao){
    x = padrao.x;
    y = padrao.y;

    return *this;
}

VETOR VETOR::operator= (double * padrao){
    x = padrao[0];
    y = padrao[1];

    return *this;
}

VETOR VETOR::operator+ (VETOR padrao){
    VETOR tmp(padrao);

    tmp.x += x;
    tmp.y += y;

    return tmp;
}

VETOR VETOR::operator+ (double* padrao){
    VETOR tmp(padrao);

    tmp.x += x;
    tmp.y += y;

    return tmp;
}

VETOR VETOR::operator- (VETOR padrao){
    VETOR tmp;

    tmp = *this;

    tmp.x -= padrao.x;
    tmp.y -= padrao.y;

    return tmp;
}

VETOR VETOR::operator- (double* padrao){
    VETOR tmp;

    tmp = *this;

    tmp.x -= padrao[0];
    tmp.y -= padrao[1];

    return tmp;
}

VETOR VETOR::operator+= (VETOR padrao){
    *this = *this + padrao;

    return *this;
}

VETOR VETOR::operator+= (double* padrao){
    *this = *this + padrao;

    return *this;
}

double VETOR::norma(){
    return sqrt(pow(x, 2) +pow(y, 2));
}

LINEAR::LINEAR(){
    Posicao = {0.0, 0.0};
    Velocidade = {0.0, 0.0};
    Aceleracao = {0.0, 0.0};
}

LINEAR::LINEAR(VETOR p, VETOR v, VETOR a){
    Posicao = p;
    Velocidade = v;
    Aceleracao = a;
}

ANGULAR::ANGULAR(){
    Posicao = {0.0, 0.0};
    Velocidade = {0.0, 0.0};
    Aceleracao = {0.0, 0.0};
}

ANGULAR::ANGULAR(VETOR p, VETOR v, VETOR a){
    Posicao = p;
    Velocidade = v;
    Aceleracao = a;
}

VETOR::~VETOR(){
}

ANGULAR::~ANGULAR(){
}

LINEAR::~LINEAR(){
}

GEOMETRIA::~GEOMETRIA(){
}

GRAO::~GRAO(){
    if (contatos != NULL)
        delete [] contatos;
    if (vizinhos != NULL)
        delete [] vizinhos;
    if (atrito != NULL)
        delete [] atrito;
    if (Reacao_Tangencial.forca != NULL)
        delete [] Reacao_Tangencial.forca;
    if (Reacao_Tangencial.contato != NULL)
        delete [] Reacao_Tangencial.contato;
}
