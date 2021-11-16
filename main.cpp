#include <iostream>
#include <random>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

static mt19937_64 gen{
    static_cast<unsigned long long>(chrono::system_clock::now().time_since_epoch().count())
};


const int N = 64;

const double sigma = 1;
const double epsilon = 1;
const double mas = 1;

const int skolko_vivodid = 500;
const int kolvo_oper = 10; //количество операций между выводами

const double Size = 5;
const double R_cut = 3; //максимальная, дистанция взаимодействия частиц

const double V_koef = 5; //максимальная скорость
const double dt = 0.001;

const int max_N_for_Maxwell = 50;

const int kolvo_write_sred_r_2 = 10000;

const bool flag_vzaimodeistvie_throw_wall = 1;

const int Skolko_discret_start_pos = pow(N, 0.34) + 1;
const int Skolko_discret_start_speed = 10;

//-----------------------
const double Start_min_dist = Size / Skolko_discret_start_pos;
const double Start_min_speed = V_koef / Skolko_discret_start_speed;

const double sigma_6 = pow(sigma, 6);


double *massiv_for_V_max = new double[skolko_vivodid];


//-----------------------
const int N_for_radial_function = 100;
const double R_for_radial_function = 4;
const double delta_R_for_radial_function  = R_for_radial_function / N_for_radial_function;


const double efectiv_R_cut = max(R_cut, R_for_radial_function); //! нужно только для подсчета следующей формулы
const int N_for_R_cut = (Size * int(efectiv_R_cut/Size) == efectiv_R_cut) ? int(efectiv_R_cut/Size) : int(efectiv_R_cut/Size) + 1;

//!--------------------------[STRUCT DOT]--------------------------------
struct dot{
    double x[3];
    double v[3];
    double f[3];
};

void create_dot(dot* a, double x_, double y, double z, double vx, double vy, double vz){
    (a->x)[0] = x_; (a->x)[1] = y; (a->x)[2] = z;
    (a->v)[0] = vx; (a->v)[1] = vy; (a->v)[2] = vz;
    for(auto &p : (a->f)) p = 0;
}

void create_dot(dot* a, double* pos, double* speed){
    for(int i = 0; i < 3; i++){
        (a->x)[i] = pos[i];
        (a->v)[i] = speed[i];
        (a->f)[i] = 0;
    }
}

ostream& operator << (ostream& out, const dot& a){
    for(const auto &p : (a.x)) out << p << " ";
    return out;
}

bool operator == (const dot& a, const dot& b){
    return ( ((a.x)[0] == (b.x)[0]) && ((a.x)[1] == (b.x)[1]) && ((a.x)[2] == (b.x)[2]) );
}



//!--------------------------[SOMETHING]----------------------------------------
double distance_pow_2(const double* const p1, const double* const p2){
    return (p1[0] - p2[0])*(p1[0] - p2[0]) +
           (p1[1] - p2[1])*(p1[1] - p2[1]) +
           (p1[2] - p2[2])*(p1[2] - p2[2]);
}

void zero_f(dot& a){
    for(auto &F : (a.f)) F = 0;
}

void update_v(dot& a){
    for(int i = 0; i < 3; i++) (a.v)[i] += (a.f)[i]/mas * dt / 2;
}



//!--------------------------[UPDATE COORDINATES]----------------------------------
void very_bad(dot& a){
    double koef = (a.v[0]*a.v[0] + a.v[1]*a.v[1] + a.v[2]*a.v[2]) / V_koef / V_koef;
    koef = pow(koef, 0.5);
    throw runtime_error("Very big speed: " + to_string(koef) + " V_koef");
    cout << "Very big speed: " + to_string(koef) + " V_koef" << endl;
}

void zerkalnie_stenki(dot& a){
    for(int i = 0; i < 3; i++){
        auto& pos = (a.x)[i];
        auto& speed = (a.v)[i];

        if(pos < -Size || pos > 2*Size) very_bad(a);

        else if(pos > Size) {
            pos = 2*Size - pos;
            speed *= -1;
        }

        else if(pos < 0) {
            pos *= -1;
            speed *= -1;
        }
    }
}

void pereodichnie_stenka(dot& a){
    for(auto& pos : a.x){
        if(pos < -Size || pos > 2*Size) very_bad(a);

        else if(pos > Size) pos -= Size;

        else if(pos < 0) pos += Size;
    }
}

void update_coord(dot& a){
    for(int i = 0; i < 3; i++) (a.x)[i] += (a.v)[i] * dt;


    if(flag_vzaimodeistvie_throw_wall) pereodichnie_stenka(a);
        else zerkalnie_stenki(a);
}

void update_coord_with_update_delta_r(dot& a, double* pos_){
    for(int i = 0; i < 3; i++) pos_[i] += (a.v)[i] * dt;

    update_coord(a);
}


//!--------------------------[START CONDITION]-----------------------------------------
void random_position(double* m, uniform_int_distribution<int>& func){
    for(int i = 0; i < 3; i++) m[i] = Size / (Skolko_discret_start_pos + 1) * func(gen);
}

void random_coord(double& m, uniform_int_distribution<int>& func){
    m = Size / (Skolko_discret_start_pos + 1) * func(gen);
}

void random_speed(double* v, uniform_int_distribution<int>& func, uniform_int_distribution<int>& sign){
    for(int i = 0; i < 3; i++) v[i] = (2*sign(gen) - 1) * Start_min_speed * func(gen);
}

void start_condition(dot* m){
    uniform_int_distribution<int> dstr_x(1, Skolko_discret_start_pos);
    uniform_int_distribution<int> dstr_v(0, Skolko_discret_start_speed);
    uniform_int_distribution<int> sign(0, 1);

    double* pos = new double[3];
    double* speed = new double[3];

    int process[9] = {int(N*0.1), int(N*0.2), int(N*0.3), int(N*0.4), int(N*0.5), int(N*0.6), int(N*0.7), int(N*0.8), int(N*0.9)}; //чисто для красоты

    for(int i = 0; i < N; i++)
    {
        for(int r = 0; r < 9; r++) if(i == process[r]) {cout << (r+1)*10 << "% "; break;}

        random_position(pos, dstr_x);
        random_speed(speed, dstr_v, sign);

        create_dot(m + i, pos, speed);

        bool is_valid = 0;
        while(!is_valid)
        {
            is_valid = 1;
            for(int j = 0; j < i; j++)
                if(m[i] == m[j]){
                    is_valid = 0;

                    for(int t = 0; t < 3; t++)
                        if((m[i].x)[t] == (m[i].x)[t]) random_coord((m[i].x)[t], dstr_x);

                    break;
                }

        }
    }

    delete[] pos;
    delete[] speed;

    cout << "end of create" << endl;
}

void create_tetraidr(dot* m){
    if(N != 4) throw invalid_argument("You create the tetraidr, but:  N != 4");

    double a = Size/2;

    double b = pow(2, 1.0/6) / sqrt(2)/2;


    create_dot(m,          b + a, -b + a, -b + a,   -V_koef, V_koef, V_koef);

    create_dot(m + 1,     -b + a,  b + a, -b + a,    V_koef, -V_koef, V_koef);

    create_dot(m + 2,      a - b, -b + a,  b + a,    V_koef, V_koef, -V_koef);

    create_dot(m + 3,      b + a, b + a, b + a,     -V_koef,  -V_koef, -V_koef);
}

void Save_condition(dot* m) {
    ofstream out;
    char filename[25];
    sprintf(filename, "Data/condition.txt");
    out.open(filename, ios::app);


    for(int i = 0; i < N; i++) {
        out << m[i] << " ";

        for(auto V : (m[i]).v) out << V << " ";
        for(auto F : (m[i]).f) out << F << " ";
        out << endl;
    }

    out << endl;

    out.close();
}

void grid_condition(dot* m) {
    uniform_int_distribution<int> dstr_v(0, Skolko_discret_start_speed);
    uniform_int_distribution<int> sign(0, 1);

    if(N != 8 && N != 27 && N != 64 && N != 125 && N != 216 && N != 343 && N != 512 && N != 729 && N != 1000)
        throw invalid_argument(to_string(N) + " - no number in third pow (or N = 1)");

    int step = 1;

    while(step*step*step != N) step++;

    double delta = Size / step;

    cout << delta << endl;

    int* pos = new int[3];
    double* coord = new double[3];
    double* speed = new double[3];

    int l = 0;

    for(pos[0] = 0; pos[0] < step; pos[0]++)
        for(pos[1]= 0; pos[1] < step; pos[1]++)
            for(pos[2] = 0; pos[2] < step; pos[2]++){
                for(int i = 0; i < 3; i++) coord[i] = (pos[i] + 0.5) * delta;

                random_speed(speed, dstr_v, sign);

                cout << "l : " << l << endl;

                create_dot(m + l, coord, speed);
                l++;
            }



    delete[] coord;
    delete[] pos;
    delete[] speed;
}

void triangle_condition(dot* m){
    uniform_int_distribution<int> dstr_v(0, Skolko_discret_start_speed);
    uniform_int_distribution<int> sign(0, 1);


    if(N != 8 && N != 64 && N != 216 && N != 512 && N != 1000)
        throw invalid_argument(to_string(N) + " - no even number in third pow");

    int step = 1;

    while(step*step*step != N) step++;

    double delta = Size / step;

    cout << delta << endl;

    int* pos = new int[3];
    double* coord = new double[3];
    double* speed = new double[3];

    int l = 0;

    for(pos[0] = 0; pos[0] < step; pos[0]++)
        for(pos[1]= 0; pos[1] < step; pos[1]++)
            for(pos[2] = 0; pos[2] < step; pos[2]++){

                coord[0] = (pos[0] + 0.5) * delta;

                if(pos[0]%2 == 0){
                    coord[1] = (pos[1] + 0.25) * delta;
                    coord[2] = (pos[2] + 0.25) * delta;
                }
                else{
                    coord[1] = (pos[1] + 0.75) * delta;
                    coord[2] = (pos[2] + 0.75) * delta;
                }

                random_speed(speed, dstr_v, sign);

                cout << "l : " << l << endl;

                create_dot(m + l, coord, speed);
                l++;
            }

}



void print_Energy(const double& a, const double& b);
void print_Impuls(double* imp);
void print_Maxwell(double V_max, unsigned int* maxw, int t);
void print_Radial_function(unsigned int* kolvo_rad_f, int t);


enum class What{
    Force,
    ForceAndEnergy,
    ForceAndMaxwell,
    ForceAndEnergyAndMaxwell,
    ForceAndEnergyAndMaxwellAndRadial
};

//!--------------------------[CALCULATE ENERGY]-------------------------------------------------
double energy_k(dot* const m){
    double f = 0;

    for(int i = 0; i < N; i++)
        for(const auto& V_ : (m[i]).v) f += mas*V_*V_/2;

    return f;
}

//!--------------------------[RECALCULATE All PAIRS]-----------------------------------------------
class CalculatePair{
public:
    CalculatePair(dot* masiv) : m(masiv), a(masiv), b(masiv) {

        R_pow_2_radial_function[0] = 0;
        for(int i = 0; i < N_for_radial_function; i++){
            Radial_function[i] = 0;
            R_pow_2_radial_function[i+1] = pow((i*delta_R_for_radial_function), 2);
        }

    }


    void do_() {
        if(flag_vzaimodeistvie_throw_wall == 0 && (sdvig[0] != 0 || sdvig[1] != 0 || sdvig[2] != 0))
            throw runtime_error("Flag for sdvig is off, but sdvig != 0");

        for(int i = 0; i < 3; i++) fictiv_pos[i] = (b->x)[i] + Size * sdvig[i];

        double r_2 = distance_pow_2((a->x), fictiv_pos);


        if(MakeRadial || r_2 <= R_pow_2_radial_function[N_for_radial_function])
            for(int i = 0; i < N_for_radial_function; i++){
                if(R_pow_2_radial_function[i] < r_2 && r_2 <= R_pow_2_radial_function[i+1]){
                    Radial_function[i] += 2;
                    break;
                }
            }



        if(r_2 > R_cut * R_cut) return;


        double r_8 = pow(r_2, 4);
        double r_14 = pow(r_2, 7);

        double f_devide_r = epsilon * (48 * sigma_6 * sigma_6 / r_14 - 24 * sigma_6 / r_8);

        for(int i = 0; i < 3; i++){
            (a->f)[i] += f_devide_r * ((a->x)[i] - fictiv_pos[i]);
            (b->f)[i] += f_devide_r * (fictiv_pos[i] - (a->x)[i]);
        }

        if(CountEnergy){
            double r_6 = pow(r_2, 3);
            double r_12 = pow(r_6, 2);

            PotentialEnergy += 4*epsilon*(sigma_6 * sigma_6 / r_12 - sigma_6 / r_6);
        }
    }

    void do_with_R_cut(){
        if(!flag_vzaimodeistvie_throw_wall){
            sdvig[0] = 0; sdvig[1] = 0; sdvig[2] = 0;
            do_();
            return;
        }

        for(sdvig[0] = - N_for_R_cut; sdvig[0] <= N_for_R_cut; sdvig[0]++)
            for(sdvig[1] = - N_for_R_cut; sdvig[1] <= N_for_R_cut; sdvig[1]++)
                for(sdvig[2] = - N_for_R_cut; sdvig[2] <= N_for_R_cut; sdvig[2]++)
                    do_();
    }

    void recalculate_(){
        const double delta_V_for_Maxwell = MakeMaxwell ? (pow(max_V_2, 0.5) / max_N_for_Maxwell) : 0;


        for(int i = 0; i < N; i++){
            zero_f(m[i]);

            if(CountEnergy) for(int j = 0; j < 3; j++)
                                impuls[j] += (m[i].x)[j];

            if(MakeMaxwell){
                double* V = m[i].v;
                double module_V = pow(V[0]*V[0]+V[1]*V[1]+V[2]*V[2], 0.5);

                int position = module_V/delta_V_for_Maxwell;
                if(position >= max_N_for_Maxwell) Maxwell[max_N_for_Maxwell-1]++;
                    else Maxwell[position]++;
            }
        }

        for(int i = 0; i < N - 1; i++)
            for(int j = i+1; j < N; j++){
                a = m + i;
                b = m + j;

                do_with_R_cut();
            }
    }


    void recalculate_(What t){
        if(t == What::ForceAndEnergyAndMaxwellAndRadial){
            MakeRadial = 1;
            skolo_bilo_Radial++;

            for(auto& n : Radial_function) n = 0;
        }


        if(t == What::ForceAndEnergy || t == What::ForceAndEnergyAndMaxwell || t == What::ForceAndEnergyAndMaxwellAndRadial) {
            CountEnergy = 1;
            PotentialEnergy = 0;
            impuls[0] = 0; impuls[1] = 0; impuls[2] = 0;
        }

        if(t == What::ForceAndMaxwell || t == What::ForceAndEnergyAndMaxwell || t == What::ForceAndEnergyAndMaxwellAndRadial){
            MakeMaxwell = 1;
            skolo_bilo_Maxwell++;
            max_V_2 = 0;

            for(auto& n : Maxwell) n = 0;

            for(int i = 0; i < N; i++){
                double* V = m[i].v;

                if(V[0]*V[0]+V[1]*V[1]+V[2]*V[2] > max_V_2) max_V_2 = V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
            }

            massiv_for_V_max[skolo_bilo_Maxwell-1] = max_V_2;
        }

        recalculate_();


        if(CountEnergy) {
            CountEnergy = 0;
            print_Energy(energy_k(m), PotentialEnergy);
            print_Impuls(impuls);
        }

        if(MakeMaxwell){
            MakeMaxwell = 0;

            print_Maxwell(max_V_2, Maxwell, skolo_bilo_Maxwell);
        }

        if(MakeRadial){
            MakeRadial = 0;

            print_Radial_function(Radial_function, skolo_bilo_Radial);
        }

    }


    bool CountEnergy = 0;
    double PotentialEnergy;
    dot* const m;
    double sdvig[3];
    double fictiv_pos[3];
    double impuls[3];
    dot* a;
    dot* b;

    //--------[ДЛЯ МАКСВЕЛЛА]-------------------------------------
    bool MakeMaxwell = 0;
    int skolo_bilo_Maxwell = 0;
    unsigned int Maxwell[max_N_for_Maxwell];
    double max_V_2;

    //--------[ДЛЯ РАДИАЛЬНОЙ ФУНКЦИИ РАСПРЕДЕЛЕНИЯ]--------------
    bool MakeRadial = 0;
    int skolo_bilo_Radial = 0;
    unsigned int Radial_function[N_for_radial_function];
    double R_pow_2_radial_function[N_for_radial_function + 1]; //границы для квадратов радиусов
};




//!--------------------------[PRINT FILE]----------------------------------------------
void print(dot* m, double t){
    ofstream out;
    char filename[25];
    sprintf(filename, "Data/full.xyz");
    out.open(filename, ios::app);

    out << N + 2 << endl << endl;

    int color_r[15] = {204,  51,  51,   0,  51, 153,   0,   0, 102, 204, 255, 255, 255, 255, 255};
    int color_b[15] = {204, 102, 255, 153, 255, 255, 204, 102,  51,   0,  51, 153, 153, 204,   0};
    int color_g[15] = { 51,   0,   0, 153, 153, 255, 255, 255, 204, 255, 204, 204, 102,  51,  51};


    for(int i = 0; i < N; i++) out << m[i] << " " <<  endl; //(100*i)%256 << " " << (1000*i)%256 << " " << (10000*i)%256 << endl;

    out << Size << " " << Size << " " << Size << endl; //" " << N << " " << N << " " << N  << endl;
    out << "0 0 0 " << endl; // N << " " << N << " " << N <<endl;

    out.close();
}

void clear_file(const string& filename){
    ofstream out;
    out.open(filename);
    out.close();
}

void print_Energy(const double& a, const double& b) {
    ofstream out;
    char filename[25];
    sprintf(filename, "Data/energy.txt");
    out.open(filename, ios::app);

    out << a << endl << b << endl;

    out.close();
}

void print_Impuls(double* imp){
    ofstream out;
    char filename[25];
    sprintf(filename, "Data/impuls.txt");
    out.open(filename, ios::app);

    out << imp[0]*imp[0] + imp[1]*imp[1] + imp[2]*imp[2] << endl;

    out.close();
}

void print_Maxwell(double V_max, unsigned int* maxw, int t){
    ofstream out;
    string buffer = "Data/Maxwell/Maxwell_" + to_string(t) + ".txt";
    out.open(buffer);

    out << sqrt(V_max) << endl;

    for(int i = 0; i < max_N_for_Maxwell; i++)
        out << maxw[i] << endl;

    out.close();
}

void print_Sred_r_2(const double& sr){
    ofstream out;
    char filename[25];
    sprintf(filename, "Data/Sred_r_2.txt");
    out.open(filename, ios::app);

    out << sr << endl;

    out.close();
}

void print_AKFC(const double beg[][3], const dot* all_dot){
    ofstream out;
    char filename[25];
    sprintf(filename, "Data/AKFC.txt");
    out.open(filename, ios::app);

    double itog = 0;

    for(int i = 0; i < N; i++)
        for(int j = 0; j < 3; j++)
            itog += beg[i][j] * ((all_dot[i]).v)[j];

    out << itog / N << endl;

    out.close();
}

void print_Radial_function(unsigned int* kolvo_rad_f, int t){
    ofstream out;
    string buffer = "Data/Radial_function/Radial_" + to_string(t) + ".txt";
    out.open(buffer);

    out << R_for_radial_function << endl;

    for(int i = 0; i < N_for_radial_function; i++)
        out << kolvo_rad_f[i] << endl;

    out.close();
}



int plotnost(dot* m, double* coord, double r){
    int n = 0;
    r = r*r;

    for(int i = 0; i < N; i++)
        if(distance_pow_2(coord, m[i].x) < r) n++;

    return n;
}



int main()
{
    cout << "N_for_R_cut " << N_for_R_cut << endl;
    cout << "Skolko_discret_start_pos " << Skolko_discret_start_pos << endl;
    clear_file("Data/full.xyz");
    clear_file("Data/energy.txt");
    clear_file("Data/condition.txt");
    clear_file("Data/total_energy.txt");
    clear_file("Data/impuls.txt");
    clear_file("Data/Sred_r_2.txt");
    clear_file("Data/AKFC.txt");




    dot* old = new dot[N];


    triangle_condition(old);
    //if(N == 4) create_tetraidr(old);
      //  else start_condition(old);


    CalculatePair updater(old);


    updater.recalculate_();
    print(old, 0);


    int K = skolko_vivodid*kolvo_oper;
    double t = 0;

    bool check_sred_r_2 = 0;
    int how_long_write_delta_r = kolvo_write_sred_r_2;
    double delta_pos[N][3];
    for(int i = 0; i < N; i++)
        for(int j = 0; j < 3; j++)
            delta_pos[i][j] = 0;

    double begin_V[N][3]; //! ЗАПИСЬ АКФС ИДЕТ ПАРАЛЕЛЬНО С Sred_r_2


    int firssssst = 0;
    while(K--)
    {
        if(firssssst){
            double proverka_V_max = 0;
            for(int i = 0; i < N; i++){
                double* ykaz_V = old[i].v;

                if(proverka_V_max < ykaz_V[0]*ykaz_V[0] + ykaz_V[1]*ykaz_V[1] + ykaz_V[2]*ykaz_V[2])
                    proverka_V_max = ykaz_V[0]*ykaz_V[0] + ykaz_V[1]*ykaz_V[1] + ykaz_V[2]*ykaz_V[2];
            }
            cout << sqrt(proverka_V_max) << endl;
            firssssst--;
        }



        //updater.recalculate_(); //!добавил отсебятину, можно удалить эту строчку безболезнено


        int step = (K%kolvo_oper == 0) ? skolko_vivodid - K/kolvo_oper : 0; //step пробегает от 1 до skolko_vivodid

        if(step == 300) {
            check_sred_r_2 = 1;

            for(int i = 0; i < N; i++)
                for(int j = 0; j < 3; j++)
                    begin_V[i][j] = (old[i].v)[j];
        }

        double sred_r_2 = 0;
        for(int i = 0; i < N; i++) {
            update_v(old[i]);

            if(check_sred_r_2){
                update_coord_with_update_delta_r(old[i], delta_pos[i]);

                sred_r_2 += delta_pos[i][0]*delta_pos[i][0] + delta_pos[i][1]*delta_pos[i][1] + delta_pos[i][2]*delta_pos[i][2];
            }
            else update_coord(old[i]);
        }

        if(check_sred_r_2 && how_long_write_delta_r--) {
            print_Sred_r_2(sred_r_2/N);
            print_AKFC(begin_V, old);
        }

        if(how_long_write_delta_r == 0) check_sred_r_2 = 0;

        if(step) updater.recalculate_(What::ForceAndEnergyAndMaxwellAndRadial);
            else updater.recalculate_();

        for(int i = 0; i < N; i++) update_v(old[i]);

        t += dt;
        if(step)
        {
            print(old, t);
            cout << step << endl;
            //Save_condition(old);
        }



    }

    //for(int i = 0; i < skolko_vivodid; i++) cout << massiv_for_V_max[i] << endl;


    delete[] old;
    //delete[] now;

    return 0;
}
