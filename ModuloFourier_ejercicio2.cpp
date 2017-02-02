// g++ -O2 ModuloFourier_ejercicio2.cpp -o test -lfftw3 -lrt `pkg-config --cflags opencv` `pkg-config --libs opencv`

// preprocesor directives
#include <opencv2/core/core.hpp>          // OpenCV
#include <opencv2/highgui/highgui.hpp>
#include <iostream>                       // cout
#include <sys/time.h>                     // funciones de tiempo
#include <math.h>                         // funciones matematicas
#include <fftw3.h>                        // lfunciones del FFTW
using namespace std;
using namespace cv;

int main(int argc, char const *argv[])
{

  //espera ruta para leer imagen
  if(argc != 2)
  {
    cout << "Usage: ./ejecutable lugarDondeEstaLaImagen" << endl;
    return -1;
  }

  //lee imagen, convierte  a tonos de gris, despliega imagen
  Mat imageIn = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);  //read file
  if(!imageIn.data)                                             //Check for invalid input
  {
    cout << "Could not open or find the image" << endl;;
    return -1;
  }
  Mat imageIn1;

  imageIn.copyTo(imageIn1); //Otra instancia de la imagenIn, aqui se modificara para imagen de salida1

  namedWindow("Imagen de entrada", WINDOW_AUTOSIZE); //Create a window for display
  imshow("Imagen de entrada",imageIn);               //Show our image inside it

  //marca tiempo de inicio
  struct timeval start ,end; //variables de tiempo
  gettimeofday(&start, NULL);

  //Calculo el tamaño de la imagen si es potencia de 2
  //si no es aumento el tamaño a una potencia
  double n = log10( double(imageIn.rows) )/log10(2.);
  double m = log10( double(imageIn.cols) )/log10(2.);

  if(double( int(n) ) != n)
    n = int(n) + 1;
  if(double( int(n) ) != m)
    m = int(m) + 1;

  int rowup = int( pow(2.,n) );
  int colup = int( pow(2.,m) );

  //separo memoria usando mallocate memory, update variables
  fftw_complex *dato0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * rowup * colup);
  fftw_complex *dato1 = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * rowup * colup );

  //copia los valores de la imagen
  //ojo paso de una imagen de (ren,col) a un vector de rowup x colup
  //la imagen son puros valores reales, por lo que imaginario es cero

  for( int r = 0; r < imageIn.rows; r++)
    for(int c = 0; c < imageIn.cols; c++)
    {
      dato0[r*colup+c][0] = double (imageIn.at<unsigned char>(r,c))*pow(-1.0,double(r+c));
      dato0[r*colup+c][1] = 0.;

      dato1[r*colup+c][0] = double (imageIn.at<unsigned char>(r,c))*pow(-1.0,double(r+c));
      dato1[r*colup+c][1] = 0.;

    }

  //calcula la transformada de Fourier
  fftw_plan process = fftw_plan_dft_2d(rowup,colup,dato0,dato0, FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute( process);
  fftw_destroy_plan( process );

  fftw_plan process1 = fftw_plan_dft_2d(rowup,colup,dato1,dato1, FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute( process1);
  fftw_destroy_plan( process1 );

  //Calcula el espectro de potencia
  Mat dummy( rowup, colup, CV_64F);
  Mat dummy1( rowup, colup, CV_64F);

  //Calcula el centro de la imagen
  int uo = rowup/2, vo = colup/2;
  double Do = 0.1, W = 0.05;

  for( int u = 0; u<rowup; u++ )
    for( int v = 0; v<colup; v++ )
    {

      //extract values
      double real = dato0[u*colup+v][0];
      double imag = dato0[u*colup+v][1];

      //Calcula la posicion con respecto al centro del dominio frecuencial
      double pu = double(u-uo)/double(rowup);
      double pv = double(v-vo)/double(colup);
      double D = sqrt(pu*pu + pv*pv);

      //filtro pasabajas
      double valor = (D < Do)?1.0:0.0; //filtro
      //double q = ( D*D - Do*Do)/(D*W);
      //double valor = 1 - exp( -q*q );
      double valor1 = valor;

      //eliminar la informacion correspondiente los componentes consenoidoles
///*

      //acotamos entre dos angulos los componentes frecuenciales deseados
      double angulo = atan2(pu,pv)*180.0/M_PI;
      if(angulo >= 75)
      {
        if(angulo <= 90)
        {
          valor *= 1.0;
        }else
        {
          valor = 0.0;
        }
      }else
      {
        if(fabs(angulo) >= 90)
        {
          if(fabs(angulo)<=105)
          {
            valor *=1.0;
          }else{
            valor = 0.0;
          }
        }else{
          valor = 0.0;
        }
      }
//*/
      dato0[u*colup + v][0] = real*valor;
      dato0[u*colup + v][1] = imag*valor;

      //Calcula el espectro de potencia
      real = dato0[u*colup+v][0];
      imag = dato0[u*colup+v][1];



      dummy.at<double>(u,v) = log10(1.0 + real * real +imag * imag);


    //extract values
      double real1 = dato1[u*colup+v][0];
    double imag1 = dato1[u*colup+v][1];

    //Calcula la posicion con respecto al centro del dominio frecuencial

    //filtro pasabajas
    //double valor1 = (D < Do)?1.0:0.0; //filtro
    //double q = ( D*D - Do*Do)/(D*W);
    //double valor = 1 - exp( -q*q );

    //acotamos entre dos angulos los contenidos frecuenciales deseados
    if(angulo >= 50)
    {
      if(angulo <= 55)
      {
        valor1 *= 1.0;
      }else
      {
        valor1 = 0.0;
      }
    }else
    {
      if(fabs(angulo)>= 125)
      {
        if(fabs(angulo)<=130)
        {
          valor1 *=1.0;
        }else{
          valor1 = 0.0;
        }
      }else{
        valor1 = 0.0;
      }
    }
    //*/
    dato1[u*colup + v][0] = real1*valor1;
    dato1[u*colup + v][1] = imag1*valor1;

    //Calcula el espectro de potencia
    real1 = dato1[u*colup+v][0];
    imag1 = dato1[u*colup+v][1];

    dummy1.at<double>(u,v) = log10(1.0 + real1 * real1 + imag1 * imag1);
}

  normalize(dummy, dummy, 0,1, CV_MINMAX);
  namedWindow("Espectro de frecuencia",WINDOW_AUTOSIZE);
  imshow("Espectro de frecuencia", dummy);

  normalize(dummy1, dummy1, 0,1, CV_MINMAX);
  namedWindow("Espectro de frecuencia1",WINDOW_AUTOSIZE);
  imshow("Espectro de frecuencia1", dummy1);

  //Calcula  la antitransformada de Fourier
  process = fftw_plan_dft_2d(rowup,colup,dato0,dato0,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute( process );
  fftw_destroy_plan( process );

  process1 = fftw_plan_dft_2d(rowup,colup,dato1,dato1,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_execute( process1 );
  fftw_destroy_plan( process1 );


  //Se obtiene la imagen filtrada
  Mat imageOut( imageIn );
  Mat imageOut1( imageIn1 );

  double cte = 1.0 / ( double(rowup)*double(colup) ); //Constante de normalizacion
  for ( int r = 0; r < imageIn.rows; r++)
    for( int c = 0; c < imageIn.cols; c++)
    {
      //extract values
      double real = cte * dato0[r*colup+c][0]*pow(-1.0,double(r+c));
      double imag = cte * dato0[r*colup+c][1]*pow(-1.0,double(r+c));

      double real1 = cte * dato1[r*colup+c][0]*pow(-1.0,double(r+c));
      double imag1 = cte * dato1[r*colup+c][1]*pow(-1.0,double(r+c));


      //guarda en variable
      double valor =sqrt(real*real+imag*imag);
      valor = valor >255.0?255.0:valor;

      double valor1 =sqrt(real1*real1 + imag1*imag1);
      valor1 = valor1 >255.0?255.0:valor1;

      imageOut.at<unsigned char>(r,c) = (unsigned char) valor;
      imageOut1.at<unsigned char>(r,c) = (unsigned char) valor1;
    }
  namedWindow("Imagen de salida",WINDOW_AUTOSIZE);
  imshow("Imagen de salida",imageOut);

  namedWindow("Imagen de salida1",WINDOW_AUTOSIZE);
  imshow("Imagen de salida1",imageOut1);

  //marca de fin de tiempo cronometrado
  gettimeofday(&end, NULL);

  //calcula el tiempo utilizado en milisegundos
  double startms = double(start.tv_sec)*1000 + double(start.tv_usec)/1000;
  double endms = double(end.tv_sec)*1000 + double(end.tv_usec)/1000;
  double ms = endms - startms;
  cout << endl << "Tiempo empleado : "<< ms << " mili-segundos" << endl;

  //termina la ejecucion del programa
  waitKey(0);  //Wait for a keystroke in the window
  fftw_free(dato0);
  return 0;
}
