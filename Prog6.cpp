#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>

using namespace std;
// Here are some standard constants and conversions
#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

const float earthradius = 3963.0;    // [miles]
const float distance_accuracy = 5.0; // [miles]

const int national_minpop = 1000000;

const float national_dist = 150.0; // [miles]
const float regional_dist = 100.0; // [miles]

const float local_maxdist = 50.0; // [miles]

const float plane_speed = 520.0; // [mph]
const float truck_speed = 65.0;  // [mph]

enum city_t
{
    LOCAL,
    METRO,
    REGIONAL,
    NATIONAL
};
enum color_t
{
    WHITE,
    BLACK
};
// Here is the city class
class city
{
private:
    // all variables that a city has
    string name;
    float latitude;
    float longitude;
    int population;
    city_t type;
    int index;

public:
    static int nameWidth; // for printing purposes
    friend istream &operator>>(istream &in, city &c)
    { // input operator to read in the .csv data
        string line;
        getline(in, line);

        // replace commas with spaces and spaces with underscores
        for (int i = 0; i < int(line.size()); i++)
        {

            if (line[i] == ',')
            {
                line[i] = ' ';
            }
            else if (line[i] == ' ')
            {
                line[i] = '_';
            }
        }

        // use stringstream to parse the line
        stringstream ss(line);
        string cityName, stateName, type;
        ss >> cityName >> stateName >> type >> c.latitude >> c.longitude >> c.population;

        // get the city type
        if (type == "LOCAL")
        {
            c.type = LOCAL;
        }
        else if (type == "METRO")
        {
            c.type = METRO;
        }
        c.latitude *= DEG2RAD;
        c.longitude *= DEG2RAD;

        // combine city and state for name
        c.name = cityName + "_" + stateName;
        nameWidth = max(nameWidth, int(c.name.size()));

        static int index = 0; // index the city
        c.index = index++;

        return in;
    }
    bool operator<(const city &c) const
    { // overload the < operator to sort them later
        return c.population < population;
    }
    friend ostream &operator<<(ostream &out, const city &c)
    { // output operator to print the city

        static unordered_map<city_t, string> type_to_string = {{LOCAL, "LOCAL"}, {METRO, "METRO"}, {REGIONAL, "REGIONAL"}, {NATIONAL, "NATIONAL"}};

        out << setw(nameWidth + 3) << left << setfill('.') << c.get_name() << setfill(' ')
            << "  " << setw(8) << left << type_to_string[c.get_type()]
            << "  " << setw(8) << right << c.get_population()
            << "  " << setw(7) << fixed << setprecision(2) << right << c.get_latitude() * RAD2DEG
            << "  " << setw(7) << fixed << setprecision(2) << right << c.get_longitude() * RAD2DEG;
        return out;
    }
    // getters and setters for the city
    string get_name() const
    {
        return name;
    }
    int get_index() const
    {
        return index;
    }
    float get_latitude() const
    {
        return latitude;
    }
    float get_longitude() const
    {
        return longitude;
    }
    int get_population() const
    {
        return population;
    }
    city_t get_type() const
    {
        return type;
    }
    void set_type(city_t t)
    {
        type = t;
    }
};

int city::nameWidth = 0; // initialize the name width to 0

// Here is the matrix class for a triangular matrix
class matrix
{
private:
    vector<float> data;
    int rows;
    int cols;

public:
    matrix(int n) : rows(n), cols(n)
    { // constructor uses the formula for a triangular matrix
        data.resize(n * (n + 1) / 2);
    }
    float &operator()(int i, int j)
    { // operator overload switches the indicies if neccessary
        if (i < j)
        {
            return data[j * (j + 1) / 2 + i];
        }
        return data[i * (i + 1) / 2 + j];
    }
};
// create vertex table reads in the .csv and indexes the cities
void create_vertex_table(char *fname, vector<city> &vertex_table)
{
    ifstream fin(fname);
    if (!fin)
    {
        cerr << "Error: could not open file " << fname << endl;
        exit(1);
    }

    city c;
    while (fin >> c) // simple bc the >> operator is overloaded
    {
        vertex_table.push_back(c);
    }
    sort(vertex_table.begin(), vertex_table.end()); // sort using the overloaded < operator
}
// update vertex table compares cities to eachother and updates their type
void update_vertex_table(vector<city> &vertex_table, matrix &dist_table)
{
    // use pointers to access the cities directly
    set<city *> national_cities;
    set<city *> regional_cities;
    for (auto &i : vertex_table)
    {
        if (i.get_type() == METRO)
        { // if the city is a metro, check if it is a national or regional city
            if (i.get_population() > national_minpop)
                i.set_type(NATIONAL);
            else
                i.set_type(REGIONAL);
        }
        if (i.get_type() == NATIONAL)
        { // if the city is national, check if it too close to another national city
            for (auto &j : national_cities)
            {
                if (dist_table(i.get_index(), j->get_index()) < national_dist)
                {
                    if (i.get_population() < j->get_population())
                    {
                        i.set_type(REGIONAL);
                    }
                    else
                    {
                        j->set_type(REGIONAL);
                        national_cities.erase(j);
                    }
                }
            }
            if (i.get_type() == NATIONAL)
                national_cities.insert(&i);
        }
        if (i.get_type() == REGIONAL)
        { // if the city is regional, check if it is too close to another regional city
            for (auto &j : regional_cities)
            {
                if (dist_table(i.get_index(), j->get_index()) < regional_dist)
                {
                    if (i.get_population() < j->get_population())
                    {
                        i.set_type(LOCAL);
                    }
                    else
                    {
                        j->set_type(LOCAL);
                        regional_cities.erase(j);
                    }
                }
            }
            if (i.get_type() == REGIONAL)
                regional_cities.insert(&i);
        }
    }
}
// create dist table uses the haversine formula to calculate the distance between two cities
void create_dist_table(vector<city> &vertex_table, matrix &dist_table)
{
    for (int i = 0; i < int(vertex_table.size()); i++)
    { // iterate through the cities and calculate the distance between them
        for (int j = i; j != int(vertex_table.size()); j++)
        {
            if (vertex_table[i].get_index() == vertex_table[j].get_index())
            {
                dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) = 0;
            }
            else
            {
                // haversine formula
                float lat1 = vertex_table[i].get_latitude();
                float lat2 = vertex_table[j].get_latitude();

                float long1 = vertex_table[i].get_longitude();
                float long2 = vertex_table[j].get_longitude();

                float dlat = (lat2 - lat1);
                float dlong = (long2 - long1);

                float a = pow(sin(dlat / 2), 2) + pow(sin(dlong / 2), 2) * cos(lat1) * cos(lat2);
                float b = 2 * asin(sqrt(a));
                float dist = earthradius * b;
                dist = distance_accuracy * round(dist / distance_accuracy);

                dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) = dist;
            }
        }
    }
}
// create time table uses the dist table to calculate the time it takes to travel between two cities
void create_time_table(vector<city> &vertex_table, matrix &dist_table, matrix &time_table)
{
    for (int i = 0; i < int(vertex_table.size()); i++)
    {
        for (int j = 0; j < int(vertex_table.size()); j++)
        {
            if (vertex_table[i].get_index() == vertex_table[j].get_index())
            { // if the cities are the same, the time is 0
                time_table(vertex_table[i].get_index(), vertex_table[j].get_index()) = 0;
            }
            else if (vertex_table[i].get_type() == NATIONAL && vertex_table[j].get_type() == NATIONAL)
            { // if the cities are national, the time is the distance divided by the plane speed
                time_table(vertex_table[i].get_index(), vertex_table[j].get_index()) = dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) / plane_speed;
            }
            else
            { // if the cities are not national, the time is the distance divided by the truck speed
                time_table(vertex_table[i].get_index(), vertex_table[j].get_index()) = dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) / truck_speed;
            }
        }
    }
}

void write_vertex_table(vector<city> &vertex_table)
{ // write the vertex table to a file
    ofstream fout("vertex_table.txt");

    for (int i = 0; i < int(vertex_table.size()); i++) // loop through and print
        fout << right << setw(4) << i << "  " << vertex_table[i] << endl;

    // close the file
    fout.close();
}
// write_dist_table writes the distance table to a file
void write_dist_table(vector<city> &vertex_table, matrix &dist_table)
{
    ofstream fout("dist_table.txt");

    vector<city> explored_cities;
    // loop through the cities and print the distance between them
    for (int i = 0; i < int(vertex_table.size()); i++)
    {
        for (int j = 0; j < i; j++)
        {
            fout << setw(4) << right << i << "  "
                 << setw(city::nameWidth + 3) << left << setfill('.') << vertex_table[i].get_name() << setfill(' ')
                 << " to "
                 << setw(city::nameWidth + 3) << left << setfill('.') << vertex_table[j].get_name() << setfill(' ')
                 << setw(8) << fixed << right << setprecision(1) << dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) << " miles" << endl;
        }
        // add the city to the explored cities
        explored_cities.push_back(vertex_table[i]);
        if (i != 0)
            fout << endl;
    }
}
// write time table writes the time table to a file
void write_time_table(vector<city> &vertex_table, matrix &time_table)
{ // pretty similar to the dist table
    ofstream fout("time_table.txt");

    vector<city> explored_cities;
    for (int i = 0; i < int(vertex_table.size()); i++)
    {
        for (int j = 0; j < i; j++)
        {
            fout << setw(4) << right << i << "  "
                 << setw(city::nameWidth + 3) << left << setfill('.') << vertex_table[i].get_name() << setfill(' ')
                 << " to "
                 << setw(city::nameWidth + 3) << left << setfill('.') << vertex_table[j].get_name() << setfill(' ')
                 << setw(8) << fixed << right << setprecision(1) << time_table(vertex_table[i].get_index(), vertex_table[j].get_index()) << " hours" << endl;
        }
        explored_cities.push_back(vertex_table[i]);
        if (i != 0)
            fout << endl;
    }
}
// write edge table writes the edge table to a file
void write_edge_table(vector<city> &vertex_table, vector<set<int>> &edge_table, matrix &dist_table, matrix &time_table)
{ // pretty similar to the dist table and time table
    ofstream fout("edge_table.txt");

    for (int i = 0; i < int(vertex_table.size()); i++)
    {
        fout << setw(4) << right << i << " " << vertex_table[i].get_name() << endl;
        for (int j : edge_table[i])
        {
            fout << setw(6) << right << j << " "
                 << setw(city::nameWidth + 3) << left << setfill('.') << vertex_table[j].get_name() << setfill(' ')
                 << setw(8) << right << fixed << setprecision(1) << dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) << " miles"
                 << setw(5) << right << fixed << setprecision(1) << time_table(vertex_table[i].get_index(), vertex_table[j].get_index()) << " hours" << endl;
        }
        if (i != int(vertex_table.size()) - 1)
            fout << endl;
    }

    fout.close();
}
// create edge table creates the edge table by linking appropriate cities
void create_edge_table(vector<city> &vertex_table, vector<set<int>> &edge_table, matrix &dist_table)
{
    // resize the edge table
    edge_table.resize(vertex_table.size());
    vector<int> local_cities, regional_cities, national_cities, non_local_cities; // create vectors to store the cities
    for (int i = 0; i < int(vertex_table.size()); i++)
    { // populate the vectors
        switch (vertex_table[i].get_type())
        {
        case LOCAL:
            local_cities.push_back(i);
            break;
        case REGIONAL:
            regional_cities.push_back(i);
            break;
        case NATIONAL:
            national_cities.push_back(i);
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < int(vertex_table.size()); i++)
    { // iterate through the cities and link them
        non_local_cities.clear();
        if (vertex_table[i].get_type() == NATIONAL)
        { // national cities are linked to all other national cities
            for (int j = 0; j < int(national_cities.size()); j++)
            {
                if (i != national_cities[j])
                    edge_table[i].insert(national_cities[j]);
            }
        }
        else if (vertex_table[i].get_type() == REGIONAL)
        { // regional cities are linked to the 3 closest non-local cities
            for (int j = 0; j < int(vertex_table.size()); j++)
            {
                if (vertex_table[j].get_type() != LOCAL && i != j)
                {
                    non_local_cities.push_back(j);
                }
            }
            // partial sort the cities to get the first 3
            partial_sort(non_local_cities.begin(), non_local_cities.begin() + 3, non_local_cities.end(), [&](int a, int b)
                         { return dist_table(vertex_table[i].get_index(), vertex_table[a].get_index()) < dist_table(vertex_table[i].get_index(), vertex_table[b].get_index()); });
            for (int j = 0; j < 3; j++)
            {
                edge_table[i].insert(non_local_cities[j]);
                edge_table[non_local_cities[j]].insert(i);
            }
        }
        else if (vertex_table[i].get_type() == LOCAL)
        { // local cities are linked to the 5 closest non-local cities and all other local cities within 50 miles
            non_local_cities.clear();
            for (int j = 0; j < int(vertex_table.size()); j++)
            {
                if (vertex_table[j].get_type() != LOCAL)
                {
                    non_local_cities.push_back(j);
                }
                else if (vertex_table[j].get_type() == LOCAL && i != j)
                {
                    if (dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) < local_maxdist)
                    {
                        edge_table[i].insert(j);
                        edge_table[j].insert(i);
                    }
                }
            }
            // partial sort the cities to get the first 5
            partial_sort(non_local_cities.begin(), non_local_cities.begin() + 5, non_local_cities.end(), [&](int a, int b)
                         { return dist_table(vertex_table[i].get_index(), vertex_table[a].get_index()) < dist_table(vertex_table[i].get_index(), vertex_table[b].get_index()); });
            for (int j = 0; j < 5; j++)
            {
                edge_table[i].insert(non_local_cities[j]);
                edge_table[non_local_cities[j]].insert(i);
            }
        }
    }
}
// here is the dijkstra route function it takes two cities and finds the shortest path between them
void dijkstra_route(string city_from, string city_to, vector<city> &vertex_table, vector<set<int>> &edge_table, int mode, matrix &dist_table, matrix &time_table, int seed = time(NULL))
{
    //declare the variables
    int index_from = -1, index_to = -1;
    if (city_from == "*" && city_to == "*")
    {// if both cities are *, pick two random cities
        srand(seed);
        city_from = vertex_table[rand() % vertex_table.size()].get_name();
        city_to = vertex_table[rand() % vertex_table.size()].get_name();
    }
    else if (city_from == "*")
    { //otherwise pick a random city for the source/destination
        srand(seed);
        city_from = vertex_table[rand() % vertex_table.size()].get_name();
    }
    if (city_to == "*")
    {
        srand(seed);
        city_to = vertex_table[rand() % vertex_table.size()].get_name();
    }
    // find the indexes of the cities using the prefix
    vector<int> index_sources, index_destinations;
    for (int i = 0; i < int(vertex_table.size()); i++)
    { //if the prefix is found, add the index to the vector
        if (vertex_table[i].get_name().substr(0, city_from.length()) == city_from)
            index_sources.push_back(i);
        if (vertex_table[i].get_name().substr(0, city_to.length()) == city_to)
            index_destinations.push_back(i);
    }
    // if the prefix are not found, print an error
    if (index_sources.size() == 0)
    {
        cout << city_from << ": prefix not found\n"
             << endl;
        return;
    }
    else if (index_destinations.size() == 0)
    {
        cout << city_to << ": prefix not found\n"
             << endl;
        return;
    }
    else
    { //otherwise, sort the prefixes to get the alphabetical first city
        sort(index_sources.begin(), index_sources.end(), [&](int a, int b)
             { return vertex_table[a].get_name() < vertex_table[b].get_name(); });
        sort(index_destinations.begin(), index_destinations.end(), [&](int a, int b)
             { return vertex_table[a].get_name() < vertex_table[b].get_name(); });

        index_from = index_sources[0];
        index_to = index_destinations[0];
    }

    // the rest is pretty much copied from the lecture notes
    vector<color_t> vColor;
    vColor.assign(vertex_table.size(), WHITE);

    vector<float> vDist;
    vDist.assign(vertex_table.size(), numeric_limits<float>::max());
    vDist[index_from] = 0;

    vector<int> vlink;
    vlink.assign(vertex_table.size(), -1);

    vlink[index_from] = index_from;

    while (1)
    {
        int next_i = -1;
        float mindist = numeric_limits<float>::max();
        for (int i = 0; i < int(vColor.size()); i++)
        {
            if (vColor[i] == WHITE && mindist > vDist[i])
            {
                next_i = i;
                mindist = vDist[i];
            }
        }
        int i = next_i;
        if (i == -1)
            break;

        vColor[i] = BLACK;

        if (i == index_to)
            break;

        for (int j : edge_table[i])
        {
            float wij = mode == 0 ? dist_table(vertex_table[i].get_index(), vertex_table[j].get_index()) : time_table(vertex_table[i].get_index(), vertex_table[j].get_index());

            if (vColor[j] == BLACK)
                continue;

            if (vDist[j] > vDist[i] + wij)
            {
                vDist[j] = vDist[i] + wij;
                vlink[j] = i;
            }
        }
    }
    // collect the cities in order to print them
    vector<int> route;
    int i = index_to;
    while (i != index_from)
    {
        route.push_back(i);
        i = vlink[i];
    }
    route.push_back(index_from);
    reverse(route.begin(), route.end()); // reverse the list to get source city first
    float total_miles = 0.00, total_hours = 0.00, current_hours = 0.00, current_miles = 0.00;
    for (int i = 0; i < int(route.size()); i++)
    {
        if (i != 0)
        { // calculate the total miles and hours
            current_hours = time_table(vertex_table[route[i]].get_index(), vertex_table[route[i - 1]].get_index());
            current_miles = dist_table(vertex_table[route[i]].get_index(), vertex_table[route[i - 1]].get_index());
            total_hours += current_hours;
            total_miles += current_miles;
        }
        cout << setw(8) << right << setprecision(2) << fixed << total_miles << " miles"
             << setw(6) << right << setprecision(2) << fixed << total_hours << " hours"
             << setw(5) << right << route[i] << " "
             << setw(city::nameWidth + 3) << left << setfill('.') << vertex_table[route[i]].get_name() << setfill(' ')
             << setw(10) << right << vertex_table[route[i]].get_population();
        if (i != 0)
        {
            cout << setw(8) << right << fixed << setprecision(2) << current_miles << " miles"
                 << setw(6) << right << fixed << setprecision(2) << current_hours << " hours";
        }
        cout << endl;
    }
    cout << endl;
}
// throw error function to print the error message
void throwError()
{
    cerr << "usage: ./Prog6 -info|dist|time [-seed=N] cities.txt" << endl;
    exit(1);
}

int main(int argc, char *argv[])
{
    // start off by cmd line parsing
    if (argc != 3 && argc != 4)
        throwError();

    bool infoMode = false;
    bool distMode = false;
    bool timeMode = false;

    string modeArgument = argv[1];
    if (modeArgument == "-info")
        infoMode = true;
    else if (modeArgument == "-dist")
        distMode = true;
    else if (modeArgument == "-time")
        timeMode = true;
    else
        throwError();

    int seed = 0;
    if (argc == 4)
    {
        string seedArgument = argv[2];
        seed = stoi(seedArgument.substr(6));
    }

    // create the vertex table
    vector<city> vertex_table;
    create_vertex_table(argv[argc - 1], vertex_table);

    // create the distance table
    matrix dist_table(vertex_table.size());
    create_dist_table(vertex_table, dist_table);

    // update the vertex table
    update_vertex_table(vertex_table, dist_table);

    // create the time table
    matrix time_table(vertex_table.size());
    create_time_table(vertex_table, dist_table, time_table);

    // create the edge table
    vector<set<int>> edge_table(vertex_table.size());
    create_edge_table(vertex_table, edge_table, dist_table);

    // if info mode is on, write the tables to files
    if (infoMode)
    {
        write_vertex_table(vertex_table);
        write_dist_table(vertex_table, dist_table);
        write_time_table(vertex_table, time_table);
        write_edge_table(vertex_table, edge_table, dist_table, time_table);
        // Stop the clock
    }

    // otherwise, run the dijkstra's loop
    else if (timeMode || distMode)
    {
        string city_from, city_to;
        while (1)
        {
            cout << "Enter> ";
            cin >> city_from >> city_to;

            if (cin.eof())
            {
                break;
            }
            if (argc == 3)
                dijkstra_route(city_from, city_to, vertex_table, edge_table, distMode ? 0 : 1, dist_table, time_table);
            else
                dijkstra_route(city_from, city_to, vertex_table, edge_table, distMode ? 0 : 1, dist_table, time_table, seed);
        }
        cout << endl;
    }
}
