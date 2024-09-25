#include <stdio.h>
#include <math.h>
#include <time.h>

//Simulação de Monte Carlo do Modelo de Ising Tridimensional

//Variáveis Globais

int a,b,c,ni,neq;// 3 Dimensões, número de iterações e número de iterações de equilibração
int seed[1];
FILE *ponteiro_video;//ponteiro global
double t,kb = 1.3806503e-23; //temperatura e constante de boltzmann
double ranf(int seed[1]);
clock_t initial_t , final_t;//conta tempo

double ranf(int seed[1])//função que gera números aleatórios
{
	    int l = 1029;
	    int c = 221591;
	    int m = 1048576;
	    seed[0] = (seed[0]*l+c)%m;
	    double fun = (1.0 * seed[0]) / (1.0*m);
	    return fun; //retorna número entre 0 e 1 em double
}


int escolhe(int p) // p é o parâmetro para escolher, pode ser a, b ou c (dimensões do sistema 3D)
{
	    int esc; //esc é o valor escolhido aleatoriamente
	    esc = round((p-1)*ranf(seed)+1);
	    return esc;
}


double energia(int vet[a+2][b+2][c+2]) //calcula a energia da configuração
{
    int i,j,k;//contadoras
    double en = 0;
    for(i=1;i<a+1;i++)
    {
        for(j=1;j<b+1;j++)
        {
            for(k=0;k<c+1;k++)//3D 
            {
                en = en -(double)(vet[i][j][k] * (vet[i-1][j][k] + vet[i+1][j][k] 
                + vet[i][j-1][k] + vet[i][j+1][k]+vet[i][j][k-1]+vet[i][j][k+1]));
            }   //acumula a energia de interação com a vizinhança
        }
    }
	
	    en = en * 0.5;
    return en;
}

void atualiza_borda (int vet[a+2][b+2][c+2]) // condicao de periodicidade
{
    int i,j;
    for (i=0;i<a+2;i++)
    {    
        for (j=0;j<b+2;j++)
        {    
            vet[i][j][0] = vet[i][j][c];//cada for e uma aresta
            vet[i][j][c+1] = vet[i][j][1];     
        }
    }
	
    for (i=0;i<a+2; i++)
    {    
        for (j=0;j<c+2;j++)
        {    
            vet[i][0][j] = vet[i][b][j];
            vet[i][b+1][j] = vet[i][1][j];
        }
    }

    for (i=0;i<b+2;i++)
    {    
        for (j=0;j<c+2;j++)    
        {    
            vet[0][i][j] =vet[a][i][j];
            vet[a+1][i][j] = vet[1][i][j];
        }
    }
    return;
}

void coleta_entrada()
{
    printf("Digite altura:\n ");
    scanf ("%d",&a);
    b = a ; 
    c = a ;
    printf("Digite a Temperatura em K:\n ");
    scanf ("%lf",&t);
    printf("Digite o numero de iteracoes\n ");
    scanf ("%d",&ni);
    return ;
}

void inicia_xadrez(int vet[a+2][b+2][c+2])
{
    int i,j,k;
    for(i=1;i<a+1;i++)//o inteiro vai de 0 até a+2, so o meio vai de 1 até a+1
    {
        //inicia os spins na configuraçao xadrez
        for(j=1;j<b+1;j++)
        {
            for(k=0;k<c+1;k++)
            {
                if((i+j+k)%2==0)
                {
                    vet[i][j][k]=1;
                }
                else
                {
                    vet[i][j][k]=-1;
                }
            }
        }
    }
    return;
}

void imprime_matriz (int vet[a+2][b+2][c+2] ,int cont)
{

    int i,j,k;
    if (cont%1==0)
    {
    fprintf (ponteiro_video,"%d\n\n",a*b*c);
        for(i=1;i<a+1;i++)    
        {
            for(j=1;j<b+1;j++)
            {
                for (k=0;k<c;k++)
                {
                    fprintf(ponteiro_video,"%d %d %d %d\n", i ,j, k , (1+vet[i][j][k])/2);
                }
            }
        }        
    }
    return;
}

double magnetizacao(int vet[a+2][b+2][c+2])
{
    int i,j,k;
    double mag=0;  
    for(i=0;i<a;i++)
    {
        for (j=0;j<b;j++)
        {
            for(k=0;k<c;k++)
            {            
            mag = mag + (double)vet[i][j][k];
            }        
        }
    }
    return mag;
}

void conta_tempo()
{
    double delta_t ;
    final_t = clock();
    printf("%lf min \n",(double)(final_t - initial_t)/(CLOCKS_PER_SEC * 60)) ;// tempo em minutos 
    return;
}

void conta(int num, int i) //número de passos e variavel contadora
{       
    if (i%num == 0 ) //imprime a cada num passos
    {
        printf(" %lf  \n",(double)100*i/ni); //retorna a porcentagem da simulação que foi executada 
        conta_tempo(); // chama a funcao que retorna o tempo em minutos
    }
    return;
}

int main() //montecarlo
{
    clock_t initial_t , final_t;
    initial_t = clock();
    FILE *ponteiro_equilibracao,*ponteiro_video,*ponteiro_temperaturas;
    ponteiro_equilibracao =fopen("equilibracao.txt","w");
    ponteiro_temperaturas = fopen("temperaturas.txt","w");
    ponteiro_video = fopen("video.xyz","w");

    double delta_t,mag_media, mag_media_quadrados,mag, suscep, media,media_quadrados,
    u,cv,p,delta_e,beta,en,nova_en,magn;
    // p é probabilidade delta e delta_e é diferenca de energia nova e velha

    a = 10;
    b = 10;
    c = 10;
    ni = 4e6;    
    //coleta_entrada();
    neq = ni * 0.5;
    //equilibracao a partir da metade da simulacao

    int i, j,aux,aux2,aux3;//spin binários,energias,nova energia
    int spin[a+2][b+2][c+2];
    
    for (beta = 0.02 ;  beta < 0.5 ; beta = beta + 0.02)
    {
        t = 1/beta;

        inicia_xadrez(spin);
        atualiza_borda(spin);
        //imprime_matriz(spin,i);

        seed[0] = 352450618;//inicializa seed
        beta = 1/(t);//sem considerar a constante de boltzman
        media = 0;
        media_quadrados=0;
        mag_media = 0;
        mag_media_quadrados = 0;
        en = energia(spin);//funcao ja implementada que calcula a energia

        for(i=0;i<ni;i++)
        {
            aux=escolhe(a);//escolhe aleatoriamente uma posiçao em a
            aux2=escolhe(b);// posiçao em b
            aux3 = escolhe(c);//posicao em c    
            spin[aux][aux2][aux3]= -spin[aux][aux2][aux3];//troca a configuraçao
            atualiza_borda(spin);
            delta_e =  -2.0*spin[aux][aux2][aux3]*(spin[aux-1][aux2][aux3]
            + spin[aux+1][aux2][aux3]+spin[aux][aux2-1][aux3]
            +spin[aux][aux2+1][aux3]+spin[aux][aux2][aux3-1]
            +spin[aux][aux2][aux3+1]);
            if(delta_e<=0)
            {
                en = energia(spin);
            }
            //se for menor nao faz nada, so muda a borda pq mudou
            else
            {
                p = exp(-beta*delta_e);
                if(ranf(seed)<=p){//sendo o num aletorio entre 0 e 1 como sendo o proprio ranf
                    en = energia(spin);
                }
                else
                {
                spin[aux][aux2][aux3]=-spin[aux][aux2][aux3];
                atualiza_borda(spin);//
                }
            }
            //conta(1e5,i);             
            if(i>neq)
            {
                media = media + en;
                media_quadrados = media_quadrados + en*en;  
                magn = magnetizacao(spin);
                mag_media = mag_media + magn;
                mag_media_quadrados = mag_media_quadrados + magn*magn;
            }

        }        

        media = media/((double)(a*b*c*neq));
        media_quadrados = media_quadrados/((double)(a*b*c*neq));
        mag_media = mag_media/((double)(a*b*c*neq));
        mag_media_quadrados= mag_media_quadrados/ ((double)(a*b*c*neq)*(double)(a*b*c*neq));

        u = media;
        cv = (media_quadrados - media*media)/(t*t);
        mag =fabs(mag_media);
        suscep = (mag_media_quadrados - mag_media * mag_media)/t;

        printf("  Beta        U       Cv         M          X \n");
        printf("%lf %lf %lf %lf %lf\n",1/t,u,cv,mag,suscep);
        fprintf(ponteiro_temperaturas,"%lf %lf %lf %lf %lf \n",1/t,u,cv,mag,suscep);

        final_t = clock();
        delta_t = (double)(final_t - initial_t) / (CLOCKS_PER_SEC *60);
        // tempo de execução em segundos
        printf("Tempo de execução: %lf min \n",delta_t);
    }
    fclose(ponteiro_equilibracao);
    fclose(ponteiro_temperaturas);
    fclose(ponteiro_video);

    return 0;
}