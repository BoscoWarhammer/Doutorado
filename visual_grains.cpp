// Allways compile linking the libraries:
// OpenGL Libraries: -lglut -lGL -lGLU
// CImg Libraries: -lpthread -lX11

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <GL/glut.h>
#include <CImg.h>

#include "pushbox.cpp"

using namespace cimg_library;

#define DEBUGPRINT(...)       printf(__VA_ARGS__);   fflush(stdout)

int j = 0, sinal = 0, mousexi = 0, mouseyi = 0, iteracoes = 0, salvar = 0;

CImgList <unsigned char> videoImg;

char msg_padrao[] = "Testando o formato de escrita!";
char tmp[200], nome[200];
char *msg = tmp;

double escala = 1.0L, minv = 1e20, maxv=-1e20, minpx=1e20, maxpx=-1e20;

PREDITOR_CORRETOR *Caixa;

int carrega_caixa(){
    Caixa = new PREDITOR_CORRETOR;

    return 0;
}

int carrega_caixa(char *name){
    FILE *Config;
    unsigned int p = 1, n, t, e = 0;
    unsigned long long int params[10];

    Config = fopen(name,"r");
    char tmpc = 0;
    do {
        tmpc = fgetc(Config);
        e += (tmpc == 32) ? 1 : 0;

    } while(tmpc != 10);
    rewind(Config);

    if (e == 1)
        fscanf(Config, "%d %d\n", &n, &t);
    else if (e >= 2){
        fscanf(Config, "%d", &n);
        fscanf(Config, "%d", &t);
        fscanf(Config, "%d", &p);
        for (int i = 0; i < e -2; i++)
            fscanf(Config, "%llu", &params[i]);
    }
    double **dados;
    dados = new double*[t];
    if (p == 1) {
        for (int i = 0; i < t; i++){
            dados[i] = new double[3];
            fscanf(Config, "%le %le %le\n", &dados[i][0], &dados[i][1], &dados[i][2]);
        }
    }
    else {
        for (int i = 0; i < t; i++){
            dados[i] = new double[10];
            fscanf(Config, "%le %le %le %le %le %le %le %le %le %le\n", &dados[i][0], &dados[i][1], &dados[i][2], &dados[i][3], &dados[i][4], &dados[i][5], &dados[i][6], &dados[i][7], &dados[i][8], &dados[i][9]);
        }
        for (int i = 0; i < n; i++){
            maxv = maxv < sqrt(dados[i][3]*dados[i][3]+dados[i][4]*dados[i][4]) ? sqrt(dados[i][3]*dados[i][3]+dados[i][4]*dados[i][4]) : maxv;
            minv = minv > sqrt(dados[i][3]*dados[i][3]+dados[i][4]*dados[i][4]) ? sqrt(dados[i][3]*dados[i][3]+dados[i][4]*dados[i][4]) : minv;
            maxpx = maxpx < dados[i][0] ? dados[i][0] : maxpx;
            minpx = minpx > dados[i][0] ? dados[i][0] : minpx;
        }
    }
    fclose(Config);

    Caixa = new PREDITOR_CORRETOR(p, n, dados);

    for (int i = 0; i < t; i++)
        delete [] dados[i];
    delete [] dados;

    return 0;
}

char * encontrar_caminho(char * name){
    unsigned int lastpos[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    char *pchar = strchr(name,'/'), caminho[200];
    while (pchar != NULL){
        for (int i = 10; i > 0; i--)
            lastpos[i] = lastpos[i-1];
        lastpos[0] = pchar-name+1;
        pchar = strchr(pchar+1, '/');
    }
    while (name[lastpos[0]+1] == ' '){
        for (int i = 1; i < 10; i++)
            lastpos[i-1] = lastpos[i];
    }
    strncpy(caminho, name, lastpos[0]);
    caminho[lastpos[0]] = '\0';

    return caminho;
}

void output(float x, float y, char *string){
    glRasterPos2f(x,y);
    for (int i = 0; string[i] != 0; i++)
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, string[i]);
}

void display(){
    long double psi, R, G, B;
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 1.0, 1.0);
    output(0, -0.5, msg);

// Imprimindo as linhas que aparecem no desenho

    glBegin(GL_LINES);
        glColor3f(1, 0, 0);
        glVertex2f(Caixa->Bordas[0].Linear.Posicao.x +Caixa->Bordas[0].Geometria.raio, Caixa->Bordas[0].Linear.Posicao.y + Caixa->Bordas[0].Geometria.raio);
        glVertex2f(Caixa->Bordas[1].Linear.Posicao.x -Caixa->Bordas[1].Geometria.raio, Caixa->Bordas[1].Linear.Posicao.y + Caixa->Bordas[1].Geometria.raio);

        glColor3f(0, 0, 1);
        glVertex2f(Caixa->Bordas[1].Linear.Posicao.x -Caixa->Bordas[1].Geometria.raio, Caixa->Bordas[1].Linear.Posicao.y + Caixa->Bordas[1].Geometria.raio);
        glVertex2d(Caixa->Bordas[3].Linear.Posicao.x -Caixa->Bordas[3].Geometria.raio, Caixa->Bordas[3].Linear.Posicao.y - Caixa->Bordas[3].Geometria.raio);

        glColor3f(0, 1, 0);
        glVertex2d(Caixa->Bordas[3].Linear.Posicao.x -Caixa->Bordas[3].Geometria.raio, Caixa->Bordas[3].Linear.Posicao.y - Caixa->Bordas[3].Geometria.raio);
        glVertex2d(Caixa->Bordas[2].Linear.Posicao.x +Caixa->Bordas[2].Geometria.raio, Caixa->Bordas[2].Linear.Posicao.y - Caixa->Bordas[2].Geometria.raio);

        glColor3f(1, 1, 0);
        glVertex2d(Caixa->Bordas[2].Linear.Posicao.x +Caixa->Bordas[2].Geometria.raio, Caixa->Bordas[2].Linear.Posicao.y - Caixa->Bordas[2].Geometria.raio);
        glVertex2d(Caixa->Bordas[0].Linear.Posicao.x +Caixa->Bordas[0].Geometria.raio, Caixa->Bordas[0].Linear.Posicao.y + Caixa->Bordas[0].Geometria.raio);
    glEnd();

    for (int k = 0; k < Caixa->numero_graos; k++){
        psi = (Caixa->Graos[k].Linear.Velocidade.norma()-minv)/(maxv-minv);
        if (isnan(psi) || isinf(psi)){
            psi = 1.0L;
        }
        if (psi < 0.5L){
            R = 2.0L*(0.5L-psi);
            G = 2.0L*(psi);
            B = 0.0L;
        } else {
            R = 0.0L;
            G = 2.0L*(1.0L-psi);
            B = 2.0L*(psi-0.5L);
        }
        if (isnan(psi) || isinf(psi)){
            R = 1.0L;
            G = 1.0L;
            B = 1.0L;
        }
        glColor3f(R,G,B);

        glBegin(GL_LINE_LOOP);
            for (int i = 0; i < 20; i++)
                glVertex2d(Caixa->Graos[k].Linear.Posicao.x + Caixa->Graos[k].Geometria.raio * cos(2 * PI * float(i) / 20), Caixa->Graos[k].Linear.Posicao.y + Caixa->Graos[k].Geometria.raio * sin(2 * PI * float(i) / 20));
        glEnd();

/*        glBegin(GL_LINE_LOOP);
            for (int i = 0; i < 9; i++)
                glVertex2d(Caixa->Graos[k].Linear.Posicao.x + 0.8 * Caixa->Graos[k].Geometria.raio * cos(Caixa->Graos[k].Angular.Posicao.x) + 0.1 * Caixa->Graos[k].Geometria.raio * cos(2 * PI * float(i) / 9), Caixa->Graos[k].Linear.Posicao.y + 0.8 * Caixa->Graos[k].Geometria.raio * sin(Caixa->Graos[k].Angular.Posicao.x) + 0.1 * Caixa->Graos[k].Geometria.raio * sin(2 * PI * float(i) / 9));
        glEnd();*/

/*        glBegin(GL_LINES);
            glVertex2d(Caixa->Graos[k].Linear.Posicao.x, Caixa->Graos[k].Linear.Posicao.y);
            glVertex2d(Caixa->Graos[k].Linear.Posicao.x +Caixa->Graos[k].Geometria.raio * cos(Caixa->Graos[k].Angular.Posicao.x), Caixa->Graos[k].Linear.Posicao.y +Caixa->Graos[k].Geometria.raio * sin(Caixa->Graos[k].Angular.Posicao.x));
        glEnd();*/
    }
    
// Plotando a condição periódica de contorno
    for (int k = 0; k < Caixa->numero_graos; k++){
        glColor3f(0.5, 0.5, 0.5);
        
        if (Caixa->Graos[k].Linear.Posicao.x > maxpx*3/4){
            glBegin(GL_LINE_LOOP);
                for (int i = 0; i < 9; i++)
                    glVertex2d(Caixa->Graos[k].Linear.Posicao.x-(maxpx-minpx) + Caixa->Graos[k].Geometria.raio * cos(2 * PI * float(i) / 9), Caixa->Graos[k].Linear.Posicao.y + Caixa->Graos[k].Geometria.raio * sin(2 * PI * float(i) / 9));
            glEnd();
        }
        if (Caixa->Graos[k].Linear.Posicao.x < minpx*3/4){
            glBegin(GL_LINE_LOOP);
                for (int i = 0; i < 9; i++)
                    glVertex2d(Caixa->Graos[k].Linear.Posicao.x+(maxpx-minpx) + Caixa->Graos[k].Geometria.raio * cos(2 * PI * float(i) / 9), Caixa->Graos[k].Linear.Posicao.y + Caixa->Graos[k].Geometria.raio * sin(2 * PI * float(i) / 9));
            glEnd();
        }
    }

    glColor3f(1.0, 1.0, 1.0);

    glBegin(GL_LINE_LOOP);
        for (int i = 0; i < 20; i++){
            glVertex2d(Caixa->Intruso.Linear.Posicao.x + Caixa->Intruso.Geometria.raio * cos(2 * PI * float(i) / 20), Caixa->Intruso.Linear.Posicao.y + Caixa->Intruso.Geometria.raio * sin(2 * PI * float(i) / 20));
        }
    glEnd();

    glBegin(GL_LINES);
        glVertex2d(Caixa->Intruso.Linear.Posicao.x, Caixa->Intruso.Linear.Posicao.y);
        glVertex2d(Caixa->Intruso.Linear.Posicao.x +Caixa->Intruso.Geometria.raio * cos(Caixa->Intruso.Angular.Posicao.x), Caixa->Intruso.Linear.Posicao.y +Caixa->Intruso.Geometria.raio * sin(Caixa->Intruso.Angular.Posicao.x));
    glEnd();

    glutSwapBuffers();
}

void redimensiona(int w, int h){
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if (w <= h){
        glOrtho(-1, 1, -(float)h/(float)w, (float)h/(float)w, -1, 1);
    }
    else{
        glOrtho(-(float)w/(float)h, (float)w/(float)h, -1, 1, -1, 1);
    }
    glMatrixMode(GL_MODELVIEW);
}

void converterOpenGLCImg(CImg<unsigned char> *outImg){
//DEBUGPRINT("Entrou no converter.\n");
    unsigned long int width = glutGet(GLUT_WINDOW_WIDTH), height = glutGet(GLUT_WINDOW_WIDTH), tamanho = width * height;
    unsigned char dataout[3*tamanho];
    GLubyte data[3*tamanho];

    glReadBuffer(GL_BACK);
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE, data);

    unsigned long int c = 0;
    for (unsigned long int i = 0; i < tamanho; i++){
        //c += i%width == 0 ? 1 : 0;
        dataout[2*tamanho +i] = data[3*i+2*c];    // Blue
        dataout[i] = data[3*i+2*c+1];               // Green
        dataout[tamanho +i] = data[3*i+2*c+2];      // Red
    }
//DEBUGPRINT("Copiou no converter.\n");
    CImg<unsigned char> tmpImg(dataout, width, height, 1, 3);
    tmpImg.mirror("y");
    *outImg = tmpImg;
/*    CImgDisplay main_disp(outImg,"Snapshot");
    while (!main_disp.is_closed() ) {
        main_disp.wait();
    }*/
//DEBUGPRINT("Saiu do converter.\n");
}

void salvar_tela(char *name){
    CImg<unsigned char> outImg;
    converterOpenGLCImg(&outImg);
    char imgname[200];
    sprintf(imgname, "%s.bmp", name);
    outImg.save(imgname);
//    sleep(1);

/*    unsigned long int width = glutGet(GLUT_WINDOW_WIDTH), height = glutGet(GLUT_WINDOW_WIDTH), tamanho = width * height;

    unsigned char *resultado, dataout[3*tamanho];
    resultado = new unsigned char [3*tamanho];
    unsigned char dataR[3*tamanho],dataG[3*tamanho],dataB[3*tamanho];
    GLubyte data[3*tamanho];

    glReadBuffer(GL_BACK);
    glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE, data);

    unsigned long int c = 0;
    for (unsigned long int i = 0; i < tamanho; i++){
        if(i%height == 0){
            c++;
        }
            dataB[2*tamanho +i+1] = data[3*i+2*c];
            dataR[i] = data[3*i+2*c+1];
            dataG[tamanho +i+1] = data[3*i+2*c+2];
            dataout[2*tamanho +i-1] = data[3*i+2*c];    // Blue
            dataout[i] = data[3*i+2*c+1];               // Green
            dataout[tamanho +i] = data[3*i+2*c+2];      // Red
//            if ((resultado[i]!=dataout[i]) || (resultado[tamanho +i]!=dataout[tamanho +i]) || (resultado[2*tamanho +i]!=dataout[2*tamanho +i])){
//DEBUGPRINT("ERROR!\n");
//            }
    }
    CImg<unsigned char> outImgR(dataR, width, height, 1, 3);
    CImg<unsigned char> outImgG(dataG, width, height, 1, 3);
    CImg<unsigned char> outImgB(dataB, width, height, 1, 3);
DEBUGPRINT("Copiando imagem.\n");
    CImg<unsigned char> outImg1(resultado, width, height, 1, 3);
DEBUGPRINT("Fez a imagem bmp.\n");
    outImgR.mirror("y");
    outImgG.mirror("y");
    outImgB.mirror("y");
    outImg.mirror("y");
    delete [] resultado;
    CImgDisplay main_dispR(outImgR,"Snapshot");
    while (!main_dispR.is_closed() ) {
        main_dispR.wait();
    }
    CImgDisplay main_dispG(outImgG,"Snapshot");
    while (!main_dispG.is_closed() ) {
        main_dispG.wait();
    }
    CImgDisplay main_dispB(outImgB,"Snapshot");
    while (!main_dispB.is_closed() ) {
        main_dispB.wait();
    }
    CImgDisplay main_disp(outImg,"Snapshot");
    while (!main_disp.is_closed() ) {
        main_disp.wait();
    }

    outImg.save("teste.bmp");
printf("salvou bmp\n");
fflush(stdout);
	outImg.save("teste.jpeg");
printf("salvou jpeg\n");
fflush(stdout);
	outImg.save("teste.png");
printf("salvou png\n");
fflush(stdout);*/
}

void salvar_tela(){
    salvar_tela("Imagem.bmp");
}

void salvar_video(){
    CImg<unsigned char> tmpImg;
    converterOpenGLCImg(&tmpImg);
    videoImg.insert(tmpImg, j-1, false);
}

void mensagem(int i){
    switch(i){
        case 0: exit(0);
        break;
        case 1: glutFullScreen();
        break;
        case 2: glScalef(2, 2, 2);
                escala /= 2.0L;
        break;
        case 3: glScalef(0.5, 0.5, 0.5);
                escala *= 2.0L;
        break;
        case 4: salvar_tela(nome);
        break;
        case 5: salvar++;
        break;
        case 10: msg = tmp;
        break;
        case 11: msg = "Imprimir teste 1";
        break;
        case 12: msg = "Imprimir teste 2";
        break;
        case 20: sinal = 0;
        break;
        case 21: sinal = 1;
        break;
        case 22: sinal = -1;
    }
    glutPostRedisplay();
}

void temporizador(int tempo){
    glutPostRedisplay();
    glutTimerFunc(50, temporizador, 0);
    if (iteracoes > 0){
        if (j < iteracoes)
            j++;
        else {
            j = 1;
//            msg = NULL;
            salvar++;
            if (salvar == 2){
                //videoImg.save_ffmpeg("video.mpeg", 50 ,2048);
                salvar++;
//                exit(0);
            }
        }
        delete Caixa;
        sprintf(tmp, "%s_simul_%010i.dat", nome,j*10000);
        carrega_caixa(tmp);
//        if (salvar == 1){
//            sprintf(tmp, "%s%i", nome,j);
//            salvar_tela(tmp);
//            salvar_video();
//        }
    }
}

void Mouse(int botao, int estado, int x, int y){
    if ((botao == GLUT_LEFT_BUTTON) && (estado == GLUT_DOWN)){
        mousexi = x;
        mouseyi = y;
    }
    if ((botao == GLUT_LEFT_BUTTON) && (estado == GLUT_UP)){
        glTranslatef(escala * float(x -mousexi) / float(glutGet(GLUT_WINDOW_WIDTH)),escala * float(mouseyi -y) / float(glutGet(GLUT_WINDOW_HEIGHT)), 0);
        mousexi = x;
        mouseyi = y;
    }
    if (estado == 4){
        glScalef(2, 2, 2);
        escala /= 2.0L;
    }
    if (estado == 3){
        glScalef(0.5, 0.5, 0.5);
        escala *= 2.0L;
    }
}

int main(int n_arg, char** args)
{
    int resultadoanimacao = 0, resultadotexto = 0, resultadofullscreen = 0, resultadotela = 0;
    if (n_arg == 2){
        sprintf(nome, "%s", args[1]);
        carrega_caixa(nome);
    }
    else if (n_arg == 3){
        sprintf(nome, "%s", args[1]);
        iteracoes = atoi(args[2]);
//        carrega_configuracao(nome);
    } else
        return 0;

    glutInit(&n_arg, args);                                     // Inicialização da função de criação de instância OpenGL pela GLUT.
    glutInitWindowSize(1024, 768);                               // Tamanho da tela inicial em pixels.
    glutInitWindowPosition(0, 0);                               // Posição inicial da tela em relação ao fundo.
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);   // Modo inicial de exibição: Máscara dupla de buffer, Matriz RGB.
    glutCreateWindow("OpenGL");                                 // Definindo o título da tela.
    glutReshapeFunc(redimensiona);                              // Redimencionando a tela e mantendo os objetos proporcionais.

    glutMouseFunc(Mouse);

    escala *= 20.0L;                                            // Setting the scale for the next visual operation
    glScalef(0.05, 0.05, 1);                                    // Rescaling to see all particles (if their diameter are near to 1
    glTranslatef(0 ,escala * -350.0L / float(glutGet(GLUT_WINDOW_HEIGHT)), 0); // Transladating to the bottom of the window
    glutDisplayFunc(display);                                   // Função de exibição da tela, rotina que propriamente desenha.
    temporizador(0);                                            // Função que atualiza a tela dependendo do tempo estipulado.

    resultadofullscreen = glutCreateMenu(mensagem);             //
    resultadotela = glutCreateMenu(mensagem);
        glutAddMenuEntry("Zoom In", 2);
        glutAddMenuEntry("Zoom Out", 3);
    resultadoanimacao = glutCreateMenu(mensagem);               // Criação do submenu de animação.
        glutAddMenuEntry("Girar para direita.", 21);            // Propriedade 1 do submenu de animação.
        glutAddMenuEntry("Girar para esquerda.", 22);           // Propriedade 2 do submenu de animação.
        glutAddMenuEntry("Parar", 20);                          // Propriedade 3 do submenu de animação.
    resultadotexto = glutCreateMenu(mensagem);                  // Criação do submenu de texto.
        glutAddMenuEntry("Teste de Menu1", 11);                 // Propriedade 1 do submenu de texto.
        glutAddMenuEntry("Teste de Menu2", 12);                 // Propriedade 2 do submenu de texto.
        glutAddMenuEntry("Padrao", 10);                         // Propriedade 3 do submenu de texto.
    glutCreateMenu(mensagem);                                   // Criação do menu principal.
        glutAddMenuEntry("FullScreen", 1);
        glutAddMenuEntry("Salvar Imagem", 4);
        glutAddMenuEntry("Salvar Video", 5);
        glutAddSubMenu("Zoom", resultadotela);
        glutAddSubMenu("Animar", resultadoanimacao);            // Propriedade 1 do menu principal.
        glutAddSubMenu("Mensagem", resultadotexto);             // Propriedade 2 do menu principal.
        glutAddMenuEntry("Sair", 0);                            // Propriedade 3 do menu principal.
    glutAttachMenu(GLUT_RIGHT_BUTTON);                          // Atribuir o evento de clicar direito no mouse com o menu principal.

    glutMainLoop();                                             // Função de loop para desenhar na tela.
    return 0;
}
