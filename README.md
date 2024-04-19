# City Routing Application

This C++ application calculates optimal routes between cities using Dijkstra's algorithm. It takes into account real-world factors such as city population, location, and transportation speeds to provide accurate routing results.

## Features

- Efficient data structures (custom matrix and city classes) to store and manipulate city data
- Parsing and processing of city data from CSV files
- Generation of distance and time tables
- Creation of an edge table representing city connections based on distance and city type (local, regional, national)
- Multiple modes of operation:
  - Information mode: outputs data tables
  - Routing mode: interactively finds shortest paths between cities
- Incorporation of real-world factors for accurate routing results

## Usage

1. Compile the source code using a C++ compiler:
   ```
   g++ -o city_routing main.cpp
   ```

2. Run the application with the desired mode and input file:
   ```
   ./city_routing -info|dist|time [-seed=N] cities.csv
   ```
   - `-info`: Information mode - outputs data tables
   - `-dist`: Routing mode - finds shortest paths based on distance
   - `-time`: Routing mode - finds shortest paths based on time
   - `-seed=N` (optional): Sets the random seed for selecting cities (default: current time)
   - `cities.csv`: CSV file containing city data

3. In routing mode, enter the source and destination cities when prompted:
   ```
   Enter> [source_city] [destination_city]
   ```
   - Use `*` as a wildcard to select random cities
   - Press `Ctrl+D` or `Ctrl+Z` (Windows) to exit the application

## File Structure

- `main.cpp`: Main source code file containing the application logic
- `cities.csv`: Example CSV file with city data
- `README.md`: This readme file

## Dependencies

- C++11 or newer
