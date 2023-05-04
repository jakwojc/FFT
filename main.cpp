#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <cmath>
#include <array>
#include <stdexcept>
#include <complex>
#include <valarray>


typedef std::complex<float> Complex;
typedef std::valarray<Complex> CArray;
const float PI = 3.141592653589793238460;

using namespace std;
ofstream file;
void fft(CArray &x)
{
	// DFT
	unsigned int N = x.size(), k = N, n;
	double thetaT = 3.14159265358979323846264338328L / N;
	Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
	while (k > 1)
	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = 1.0L;
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] += x[b];
				x[b] = t * T;
			}
			T *= phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
	//// Normalize (This section make it not working correctly)
	//Complex f = 1.0 / sqrt(N);
	//for (unsigned int i = 0; i < N; i++)
	//	x[i] *= f;
}
//klasa Packet
class Packet{
protected:
string device;
string description;
long date;
public:
Packet(string Device, string Description, long Date):
device(Device) , description(Description) , date(Date) {};
virtual void abstrakcja()=0;
virtual void toString(void){
cout<<"Device:"<<device<<endl;
cout<<"Description:"<<description<<endl;
cout<<"Date:"<<date<<endl;
};
virtual ~Packet() {};
};
//klasa Sequence
template <class T, int I>
class Sequence : public Packet{
protected:
int channelNr;
string unit;
double resolution;
T buffer[I];
public:
Sequence<T,I>(string Device, string Description, long Date, int ChannelNr, string Unit, double Resolution):
Packet(Device, Description, Date), channelNr(ChannelNr), unit(Unit), resolution(Resolution){}
virtual void abstrakcja()=0;
virtual void toString(void){
Packet::toString();
cout<<"ChannelNr:"<<channelNr<<endl;
cout<<"Unit:"<<unit<<endl;
cout<<"Resolution:"<<resolution<<endl;
};
virtual ~Sequence<T,I>() {};
//zabezpieczenie przed niepoprawnym indeksem
T& operator[](int i){

if((i < 0) || (i >= I))
throw i;
return buffer[i];
}

};
//klasa TimeHistory
template <class T, int I=4000>
class TimeHistory: public Sequence<T,I> {
protected:
double sensitivity;
T buffer[I];
float realis[I];
float immaginaris[I];
float module[I];
public:
TimeHistory<T,I>(string Device, string Description, long Date, int ChannelNr, string Unit, double Resolution, double
Sensitivity):
Sequence<T,I>(Device, Description, Date, ChannelNr, Unit, Resolution), sensitivity(Sensitivity){
};
void abstrakcja() {};
void toString(void){
Sequence<T,I>::toString();
cout<<"Sensitivity:"<<sensitivity<<endl;
cout<<endl;
}
virtual ~TimeHistory<T,I>() {};
//funkcja Gauss
void Gauss(double average, double std){
std::random_device rand{};
std::mt19937 generate{rand()};
std::normal_distribution<> dist{average,std};
try{
for (int i=0; i<I; ++i){
buffer[i]=dist(generate);
}
}
catch(int i){
cerr << "Invalid index" << i << endl;
}
catch(...){
cerr << "Unknown exception";
}
}
void Impulse(){
    for (int i = 0; i<I ; i++){
    buffer[i] = 0;
    }
    buffer[10]=100;
}














//funkcja FFT
void FFT(int Sraka){
    /*//I = length
    if (I == 1) {
        realis[0] = buffer[0];
        immaginaris[0]=0;
        module[0] = realis[0];
        return;
    }
    TimeHistory<float,I/2> X_even;
    TimeHistory<float,I/2 + (I % 2 != 0)> X_odd;
    //przypisanie do xeven i xodd wartosci parzystych i nieparzystych
    for (int i=0;i<I;i++){
        if (i%2==0){
            X_even[i/2] = buffer[i];
        } else {
            X_odd[i/2] = buffer[i];
        }
    }
    //rekursja
    X_even.FFT(x);
    X_odd.FFT(x);
    float *factor_re = new float[I], factor_im = new float[I];
    //przypisywanie wartosci do wektorow factor
    for (int i=0;i<I;i++){
        factor_re[i] = cos(2*PI*i/I);
        factor_im[i] = sin(2*PI*i/I);
    }
    // teraz cos tam sie dzieje fft sie robi nie wiem o co chodzi
    for (int i=0;i<I;i++){

    }*/
    Complex *data = new Complex[I];
    for (int i=0;i<I;i++){
        data[i].real(buffer[i]);
    }
    CArray x(data,I);
    fft(x);
    for (int i=0;i<I;i++){
        realis[i] = x[i].real();
        immaginaris[i] = x[i].imag();
        module[i] = sqrt(realis[i]*realis[i] + immaginaris[i] * immaginaris[i]);
    }

}
//zapis do pliku zewnetrznego
void writetofile(){
if(file.is_open()){
file << "Time Values:" << endl;
for(int i=0; i<I;i++){
file << buffer[i] << " ";
}
file << endl;
file << "Complex Spectrum:" << endl;
for (int i=0; i<I; i++){
file << realis[i] << " + ";
file << immaginaris[i] << " i\t";
}
file << endl;
file << "Magnitude Spectrum:" << endl;
for(int i=0; i<I;i++){
file << module[i] << " ";
}
file << endl << endl << endl;
}
else{
string wyjatek = "File not opened";
throw wyjatek;
}
}
//zaprzyjaznione RMS oraz operatory przeciazone
template<class T1, int I1>
friend T1 RMS(TimeHistory<T1,I1> inp);
template<class T1, int I1>
friend TimeHistory<T1,I1> operator+(const TimeHistory<T1,I1>& add1, const TimeHistory<T1,I1>& add2);
template<class T1, int I1>
friend TimeHistory<T1,I1> operator/(const TimeHistory<T1,I1>& divide1, const TimeHistory<T1,I1>& divide2);
template<class T1, int I1>
friend void operator>=(TimeHistory<T1,I1>& sth1, TimeHistory<T1,I1>& sth2);
};
//klasa Spectrum
template <class T, int I>
class Spectrum: public Sequence<T,I> {
protected:
string skala;
int scaling;
public:
Spectrum<T,I>(string Device, string Description, long Date, int ChannelNr, string Unit, double Resolution, int Scaling):
Sequence<T,I>(Device, Description, Date, ChannelNr, Unit, Resolution), scaling(Scaling){
};
void abstrakcja() {};
void toString(void){
ChooseScaling();
Sequence<T,I>::toString();
cout<<"Scaling:"<<skala<<endl;
cout<<endl;
}
//funkcja nadzorujaca poprawny wybor skalowania
void ChooseScaling(){
switch(scaling){
case 0:
skala = "linear";
break;
case 1:
skala = "logarithmic";
break;
default:
skala = "Choose 0 for linear scale or 1 for logarithmic scaling" ;
break;
}
}
virtual ~Spectrum<T,I>() {};
};
//klasa alarm
class Alarm : public Packet{
private:
int channelNr;
double threshold;
int direction;
string kierunek;
public:
Alarm(string Device, string Description, long Date, int ChannelNr, double Threshold, int Direction):
Packet(Device, Description, Date), channelNr(ChannelNr),threshold(Threshold),direction(Direction){
};
void abstrakcja() {};
void toString(void){
ChooseDirection();
Packet::toString();
cout<<"ChannelNr:"<<channelNr<<endl;
cout<<"Threshold:"<<threshold<<endl;
cout<<"Direction:"<<kierunek<<endl;
cout<<endl;
}
//funkcja nadzorujaca poprawny wybor kierunku
void ChooseDirection(){

switch(direction){
case -1:
kierunek = "down";
break;
case 0:
kierunek = "any";
break;
case 1:
kierunek = "up";
break;
default:
kierunek = "Choose -1 for down, 0 for any, 1 for up direction" ;
break;
}
}

};
//operator zamiany
template<class T1, int I1>
void operator>=(TimeHistory<T1,I1>& sth1, TimeHistory<T1,I1>& sth2){
T1 wymien1[I1];
T1 wymien2[I1];
for(int i=0;i<I1;i++){
wymien1[i]=sth1.buffer[i];
wymien2[i]=sth2.buffer[i];
sth2.buffer[i]=wymien2[i];
sth1.buffer[i]=wymien1[i];
}
}
//operator dodawania
template<class T1, int I1>
TimeHistory<T1,I1> operator+(const TimeHistory<T1,I1>& add1, const TimeHistory<T1,I1>& add2){
TimeHistory<T1,I1> add("Device1", "add", 441, 5, "unit", 0, 0);
for(int i=0;i<I1;i++){
add.buffer[i]=add1.buffer[i]+add2.buffer[i];
}
return add;
}
//operator dzielenia
template<class T1, int I1>
TimeHistory<T1,I1> operator/(const TimeHistory<T1,I1>& divide1, const TimeHistory<T1,I1>& divide2){
TimeHistory<T1,I1> divide("Device2", "divide", 121, 1, "unit", 0, 0);
for(int i=0;i<I1;i++){
if (divide2.buffer[i]==0){
string wyjatek = "Select a divisor other than zero";
throw wyjatek;
}
divide.buffer[i]=divide1.buffer[i]/divide2.buffer[i];
}
return divide;
}
// funkcja RMS
template<class T1, int I1>
T1 RMS(TimeHistory<T1,I1> inp){
int a = I1;
double b= 0;
double c = 0;
double sum = 0;
double average = 0;
double rms = 0;
for(int i=0;i<a;i++){
double d=inp.buffer[i];
sum=sum+d;
}
average=sum/a;
for(int i=0;i<a;i++){
double e=inp.buffer[i]-average;
b=pow(e,2.0);
c=c+b;
}
rms=sqrt((c/(a-1)));
cout<< "RMS = " <<rms <<endl;
return rms;
}












int main() {
Alarm s("Device1","Alarm s",723661,3,0.5,1);
s.toString();
TimeHistory<int, 5> x("Device2","TimeHistory x",2372137,5,"unit",1660,50);
x.toString();
Spectrum<int, 7> y("Device3","Spectrum y",897123,1,"unit",1980,1);
y.toString();
try{
//otwarcie pliku do ktorego zapisane zostana dane
file.open("results.txt");
//przydzielanie i zwalnianie pamieci
try{
//inicjalizacja kanalow
TimeHistory<float,1024> test1("Device1","Sygnal 1", 50, 1, "unit", 200, 20);
TimeHistory<float,1024> test2("Device2","Sygnal 2", 120, 4, "unit", 15, 5);
TimeHistory<float,1024> test3("Device1","Sygnal 3", 34, 3, "unit", 60, 33);
TimeHistory<float,1024> test4("Device3","Sygnal 4", 2000, 2, "unit", 8, 6);
//przypisanie danych z bialym szumem
test1.Gauss(3, 0.3);
test2.Gauss(6, 0.3);
test3.Gauss(8, 0.1);
test4.Impulse();
//operacja zamiany kanalu1 i kanalu2 przy pomocy operatora przeciazonego
test1 >= test2;
file << "Channel1 and Channel2 - swap" << endl;
file << endl;
//operacja dodania kanalu2 i kanalu3 przy pomocy operatora przeciazonego
TimeHistory<float,1024> test23= test2 + test3;
file << "Channel2 and Channel3 - addition" << endl;
test23.writetofile();
//operacja dzielenia kanalu1 i kanalu4 przy pomocy operatora przeciazonego
try
{
TimeHistory<float,1024> test14= test1 / test4;
file << "Channel1 and Channel4 - division" << endl;
test14.writetofile();
}
catch(string wyjatek)
{
cout<< "Exception:" << wyjatek << endl;
}
//wyznaczenie widma zespolonego oraz amplitudowego kazdego z kanalow
test1.FFT(1);
test2.FFT(1);
test3.FFT(1);
test4.FFT(1);
//zapis danych do pliku
try
{
file << "Signal1" << endl;
test1.writetofile();
file << "Signal2" << endl;
test2.writetofile();
file << "Signal3" << endl;
test3.writetofile();
file << "Signal4" << endl;
test4.writetofile();
}
catch(string wyjatek)
{
cout<< "Exception:" << wyjatek << endl;
file.close();
}
//wypisanie na konsole wartosci RMS kazdego z kanalow
cout<<"Signal1: " <<endl;
RMS(test1);
cout<<"Signal2: " <<endl;
RMS(test2);
cout<<"Signal3: " <<endl;
RMS(test3);
cout<<"Signal4: " <<endl;
RMS(test4);
}
catch (std::runtime_error& err)
{
cout << "Exception: " << err.what();
}
}
catch(std::bad_alloc& badall)
{
cout << "Exception: " << badall.what();
}
return 0;
}