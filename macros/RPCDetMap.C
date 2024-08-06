#include <iostream>
#include <vector>
#include <cmath>

class MuonDetector {
public:
    MuonDetector() {
        barrel_eta_max = 1.2;
        endcap_eta_min = 1.2;
        endcap_eta_max = 2.4;
    }

    void get_barrel_info(double eta) {
        if (fabs(eta) < barrel_eta_max) {
            std::vector<int> wheels = {-2, -1, 0, 1, 2};
            std::vector<std::string> stations = {"MB1", "MB2", "MB3", "MB4"};
            int sectors_per_wheel = 12;
            print_barrel_info(wheels, stations, sectors_per_wheel);
        } else {
            std::cout << "Eta fora do intervalo do barril.\n";
        }
    }

    void get_endcap_info(double eta) {
        if (fabs(eta) > endcap_eta_min && fabs(eta) < endcap_eta_max) {
            std::vector<std::string> stations = {"ME1", "ME2", "ME3", "ME4"};
            int disks_per_station = 3;  // ou 4 dependendo da estação
            int sectors_per_disk = 12;
            print_endcap_info(stations, disks_per_station, sectors_per_disk);
        } else {
            std::cout << "Eta fora do intervalo dos endcaps.\n";
        }
    }

    void get_detector_info(double eta) {
        if (fabs(eta) < barrel_eta_max) {
            get_barrel_info(eta);
        } else if (fabs(eta) > endcap_eta_min && fabs(eta) < endcap_eta_max) {
            get_endcap_info(eta);
        } else {
            std::cout << "Eta fora do intervalo de detecção.\n";
        }
    }

private:
    double barrel_eta_max;
    double endcap_eta_min;
    double endcap_eta_max;

    void print_barrel_info(std::vector<int> wheels, std::vector<std::string> stations, int sectors_per_wheel) {
        std::cout << "Wheels: ";
        for (int wheel : wheels) {
            std::cout << wheel << " ";
        }
        std::cout << "\nStations: ";
        for (const std::string& station : stations) {
            std::cout << station << " ";
        }
        std::cout << "\nSectors per wheel: " << sectors_per_wheel << "\n";
    }

    void print_endcap_info(std::vector<std::string> stations, int disks_per_station, int sectors_per_disk) {
        std::cout << "Stations: ";
        for (const std::string& station : stations) {
            std::cout << station << " ";
        }
        std::cout << "\nDisks per station: " << disks_per_station << "\nSectors per disk: " << sectors_per_disk << "\n";
    }
};

// Exemplo de uso
int main() {
    double eta = 1.5;
    MuonDetector detector;
    detector.get_detector_info(eta);
    return 0;
}
