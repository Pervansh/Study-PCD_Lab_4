#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>

// #define USE_INPUT_1
#define USE_INPUT_2
#define USE_INPUT_3

#define INPUT_FILE_ID 1

const std::string input1FileName = "C:\\Users\\Alex\\Documents\\Visual Studio Projects\\BMSTU\\Semester 7\\Programm Complex Designing\\PCD_Lab_3\\PCD_Lab_3\\input1.txt";
const std::string input2FileName = "C:\\Users\\Alex\\Documents\\Visual Studio Projects\\BMSTU\\Semester 7\\Programm Complex Designing\\PCD_Lab_3\\PCD_Lab_3\\input2.txt";
const std::string input3FileName = "input3.txt";
const std::string outputFileName = "C:\\Users\\Alex\\Documents\\Visual Studio Projects\\BMSTU\\Semester 7\\Programm Complex Designing\\PCD_Lab_3\\PCD_Lab_3\\output.txt";

struct Point {
    double x;
    double y;
    double z;

    // Далее идут линейные операции для Point
    friend Point operator+ (const Point& a, const Point& b) {
        return { a.x + b.x, a.y + b.y, a.z + b.z };
    }

    friend Point operator- (const Point& a, const Point& b) {
        return { a.x - b.x, a.y - b.y, a.z - b.z };
    }

    friend Point operator* (double k, const Point& a) {
        return { k * a.x, k * a.y, k * a.z };
    }

    friend Point operator/ (const Point& a, double d) {
        return { a.x / d, a.y / d, a.z / d };
    }
};

// Файл данных ввода
struct FileInput {
    int NP;
    std::vector<Point> ps;
    int NE1, NE2;
    int type;
};

// Данные для файла вывода
struct ResultMesh {
    int NE, NP, NC;
    std::vector<std::vector<int>> elemNodeInds;    // вектор номеров узлов для каждого элемента (нумерация узлов с 0)
    std::vector<Point> nodes;                      // вектор координат узлов
    std::vector<std::vector<int>> contourNodeInds; // вектор номеров узлов для каждого контура (нумерация узлов с 0)
};

// Поиск пересечения прямых через решение СЛАУ (без проверки копланарности и т.п.)
Point lineIntersection(Point a1, Point a2, Point b1, Point b2) {
    double a = a2.x - a1.x;
    double b = b1.x - b2.x;
    double c = a2.y - a1.y;
    double d = b1.y - b2.y;
    double f = b1.x - a1.x;
    double g = b1.y - a1.y;

    double det = a * d - b * c;

    double t0 = (f * d - b * g) / det;
    double s0 = (g * a - f * c) / det;

    return { a1.x + t0 * a, a1.y + t0 * c, a1.z + t0 * (a2.z - a1.z)};
}

std::unique_ptr<FileInput> readFile(std::istream& file) {
    auto fileInput = std::make_unique<FileInput>();
    
    file >> fileInput->NP;
    fileInput->ps.resize(fileInput->NP);

    for (auto& p : fileInput->ps) {
        file >> p.x >> p.y >> p.z;
    }

    if (fileInput->NP == 2) {
        file >> fileInput->NE1;
    } else if (fileInput->NP == 4) {
        file >> fileInput->NE1 >> fileInput->NE2;
    } else {
        std::cerr << "Entered NP is not expected. It should have value of 2 or 4!" << std::endl;
    }

    file >> fileInput->type;

    return std::move(fileInput);
}

ResultMesh linElem1dMesh(const FileInput& fileInput) {
    ResultMesh mesh;
    mesh.NE = fileInput.NE1;
    mesh.NP = mesh.NE + 1;
    mesh.NC = 2;

    Point p = fileInput.ps[0]; // Начало сетки
    Point l = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы

    mesh.nodes = std::vector<Point>(mesh.NP);
    for (int i = 0; i < mesh.NP; ++i) {
        mesh.nodes[i] = i * l + p;
    }

    mesh.elemNodeInds = std::vector<std::vector<int>>(mesh.NE);
    for (int i = 0; i < mesh.NE; ++i) {
        mesh.elemNodeInds[i] = { i, i + 1 };
    }

    mesh.contourNodeInds = { {0}, {mesh.NP - 1} };

    return mesh;
}

ResultMesh quadElem1dMesh(const FileInput& fileInput) {
    ResultMesh mesh;
    mesh.NE = fileInput.NE1;
    mesh.NP = 2 * mesh.NE + 1;
    mesh.NC = 2;

    Point p = fileInput.ps[0]; // Начало сетки
    Point l = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы

    mesh.nodes = std::vector<Point>(mesh.NP);
    for (int i = 0; i < mesh.NP; ++i) {
        mesh.nodes[i] = i * l + p;
    }

    mesh.elemNodeInds = std::vector<std::vector<int>>(mesh.NE);
    for (int i = 0; i < mesh.NE; i++) {
        int j = 2 * i; // индекс 1-го узла в элементе
        mesh.elemNodeInds[i] = { j, j + 1, j + 2 };
    }

    mesh.contourNodeInds = { {0}, {mesh.NP - 1} };

    return mesh;
}

ResultMesh tetraElem2dMesh(const FileInput& fileInput) {
    ResultMesh mesh;
    mesh.NE = fileInput.NE1 * fileInput.NE2;
    mesh.NP = (fileInput.NE1 + 1) * (fileInput.NE2 + 1);
    mesh.NC = 1;

    // Угловые точки области (сетки) (против часовой стрелки: p1 -> p2 -> p4 -> p3 -> p1)
    Point p1 = fileInput.ps[0];
    Point p2 = fileInput.ps[1];
    Point p3 = fileInput.ps[3];
    Point p4 = fileInput.ps[2];

    // В-ы, соединяющие соседние узлы, для соответствующих сторон
    Point l1 = (p3 - p1) / fileInput.NE1;
    Point l2 = (p4 - p2) / fileInput.NE1;
    Point q1 = (p2 - p1) / fileInput.NE2;
    Point q2 = (p4 - p3) / fileInput.NE2;

    //Point l1 = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы
    //Point l1 = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы

    int n = fileInput.NE2 + 1;
    int m = fileInput.NE1 + 1;

    //mesh.nodes = std::vector<Point>();
    mesh.nodes.reserve(mesh.NP);
    for (int i = 0; i <= fileInput.NE1; ++i) {
        for (int j = 0; j <= fileInput.NE2; ++j) {
            auto intersection = lineIntersection(p1 + i * l1, p2 + i * l2, p1 + j * q1, p3 + j * q2);
            mesh.nodes.push_back(intersection);
        }
    }

    //mesh.elemNodeInds = std::vector<std::vector<int>>();
    mesh.elemNodeInds.reserve(mesh.NE);
    for (int i = 0; i < fileInput.NE1; ++i) {
        for (int j = 0; j < fileInput.NE2; ++j) {
            // номера идут против часовой стрелки
            mesh.elemNodeInds.push_back({ i * n + j, i * n + j + 1, (i + 1) * n + j + 1, (i + 1) * n + j });
        }
    }

    // mesh.contourNodeInds = { {0, fileInput.NE1, mesh.NP - 1, mesh.NP - fileInput.NE1} };
    //mesh.contourNodeInds = { std::vector<int>(2 * (fileInput.NE1 + fileInput.NE2)) };
    std::vector<int> contourNodeInds;
    contourNodeInds.reserve(2 * (fileInput.NE1 + fileInput.NE2));
    for (int j = 0; j < fileInput.NE2; ++j) {
        contourNodeInds.push_back(j);
    }

    for (int i = 0; i < fileInput.NE1; ++i) {
        contourNodeInds.push_back(i * n + n - 1);
    }

    for (int j = 0; j < fileInput.NE2; ++j) {
        contourNodeInds.push_back((m - 1) * n + n - 1 - j);
    }

    for (int i = 0; i < fileInput.NE1; ++i) {
        contourNodeInds.push_back((m - 1 - i) * n);
    }

    mesh.contourNodeInds = { std::move(contourNodeInds) };

    return mesh;
}

ResultMesh diagElem2dMesh(const FileInput& fileInput, bool alterDiag = false) {
    ResultMesh mesh;
    mesh.NE = 2 * fileInput.NE1 * fileInput.NE2;
    mesh.NP = (fileInput.NE1 + 1) * (fileInput.NE2 + 1);
    mesh.NC = 1;

    // Угловые точки области (сетки)
    Point p1 = fileInput.ps[0];
    Point p2 = fileInput.ps[1];
    Point p3 = fileInput.ps[3];
    Point p4 = fileInput.ps[2];

    // В-ы, соединяющие соседние узлы, для соответствующих сторон
    Point l1 = (p3 - p1) / fileInput.NE1;
    Point l2 = (p4 - p2) / fileInput.NE1;
    Point q1 = (p2 - p1) / fileInput.NE2;
    Point q2 = (p4 - p3) / fileInput.NE2;

    //Point l1 = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы
    //Point l1 = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы

    int n = fileInput.NE2 + 1;
    int m = fileInput.NE1 + 1;

    //mesh.nodes = std::vector<Point>();
    mesh.nodes.reserve(mesh.NP);
    for (int i = 0; i <= fileInput.NE1; ++i) {
        for (int j = 0; j <= fileInput.NE2; ++j) {
            auto intersection = lineIntersection(p1 + i * l1, p2 + i * l2, p1 + j * q1, p3 + j * q2);
            mesh.nodes.push_back(intersection);
        }
    }

    //mesh.elemNodeInds = std::vector<std::vector<int>>();
    mesh.elemNodeInds.reserve(mesh.NE);
    for (int i = 0; i < fileInput.NE1; ++i) {
        int n = fileInput.NE2 + 1;
        for (int j = 0; j < fileInput.NE2; ++j) {
            // номера идут против часовой стрелки
            if (alterDiag) {
                mesh.elemNodeInds.push_back({ i * n + j, i * n + j + 1, (i + 1) * n + j });
                mesh.elemNodeInds.push_back({ i * n + j + 1, (i + 1) * n + j + 1, (i + 1) * n + j });
            } else {
                mesh.elemNodeInds.push_back({ i * n + j, (i + 1) * n + j + 1, (i + 1) * n + j });
                mesh.elemNodeInds.push_back({ i * n + j, i * n + j + 1, (i + 1) * n + j + 1 });
            }
        }
    }

    // mesh.contourNodeInds = { {0, fileInput.NE1, mesh.NP - 1, mesh.NP - fileInput.NE1} };

    std::vector<int> contourNodeInds;
    contourNodeInds.reserve(2 * (fileInput.NE1 + fileInput.NE2));
    for (int j = 0; j < fileInput.NE2; ++j) {
        contourNodeInds.push_back(j);
    }

    for (int i = 0; i < fileInput.NE1; ++i) {
        contourNodeInds.push_back(i * n + n - 1);
    }

    for (int j = 0; j < fileInput.NE2; ++j) {
        contourNodeInds.push_back((m - 1) * n + n - 1 - j);
    }

    for (int i = 0; i < fileInput.NE1; ++i) {
        contourNodeInds.push_back((m - 1 - i) * n);
    }

    mesh.contourNodeInds = { std::move(contourNodeInds) };

    return mesh;
}

ResultMesh crossElem2dMesh(const FileInput& fileInput) {
    ResultMesh mesh;
    mesh.NE = 4 * fileInput.NE1 * fileInput.NE2;
    mesh.NP = (fileInput.NE1 + 1) * (fileInput.NE2 + 1) + fileInput.NE1 * fileInput.NE2;
    mesh.NC = 1;

    // Угловые точки области (сетки)
    Point p1 = fileInput.ps[0];
    Point p2 = fileInput.ps[1];
    Point p3 = fileInput.ps[3];
    Point p4 = fileInput.ps[2];

    // В-ы, соединяющие соседние узлы, для соответствующих сторон
    Point l1 = (p3 - p1) / fileInput.NE1;
    Point l2 = (p4 - p2) / fileInput.NE1;
    Point q1 = (p2 - p1) / fileInput.NE2;
    Point q2 = (p4 - p3) / fileInput.NE2;

    //Point l1 = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы
    //Point l1 = (fileInput.ps[1] - fileInput.ps[0]) / (mesh.NP - 1); // в-р, соединяющий соседние узлы

    int n = fileInput.NE2 + 1;
    int m = fileInput.NE1 + 1;

    // Добавления угловых узлов элементов (без "центральных")
    mesh.nodes.reserve(mesh.NP);
    for (int i = 0; i <= fileInput.NE1; ++i) {
        for (int j = 0; j <= fileInput.NE2; ++j) {
            auto intersection = lineIntersection(p1 + i * l1, p2 + i * l2, p1 + j * q1, p3 + j * q2);
            mesh.nodes.push_back(intersection);
        }
    }

    // Добавление элементов и "центральных" узлов элементов
    mesh.elemNodeInds.reserve(mesh.NE);
    for (int i = 0; i < fileInput.NE1; ++i) {
        for (int j = 0; j < fileInput.NE2; ++j) {
            Point lp1 = mesh.nodes[i * n + j];
            Point lp2 = mesh.nodes[(i + 1) * n + j];
            Point lp3 = mesh.nodes[i * n + j + 1];
            Point lp4 = mesh.nodes[(i + 1) * n + j + 1];

            Point midP = lineIntersection(lp1, lp4, lp2, lp3);

            mesh.nodes.push_back(midP);
            int last = mesh.nodes.size() - 1;

            // номера идут против часовой стрелки
            mesh.elemNodeInds.push_back({ i * n + j, i * n + j + 1, last });
            mesh.elemNodeInds.push_back({ i * n + j + 1, (i + 1) * n + j + 1, last });
            mesh.elemNodeInds.push_back({ (i + 1) * n + j + 1, (i + 1) * n + j, last });
            mesh.elemNodeInds.push_back({ (i + 1) * n + j, i * n + j, last });
        }
    }

    // mesh.contourNodeInds = { {0, fileInput.NE1, mesh.NP - 1, mesh.NP - fileInput.NE1} };

    std::vector<int> contourNodeInds;
    contourNodeInds.reserve(2 * (fileInput.NE1 + fileInput.NE2));
    for (int j = 0; j < fileInput.NE2; ++j) {
        contourNodeInds.push_back(j);
    }

    for (int i = 0; i < fileInput.NE1; ++i) {
        contourNodeInds.push_back(i * n + n - 1);
    }

    for (int j = 0; j < fileInput.NE2; ++j) {
        contourNodeInds.push_back((m - 1) * n + n - 1 - j);
    }

    for (int i = 0; i < fileInput.NE1; ++i) {
        contourNodeInds.push_back((m - 1 - i) * n);
    }

    mesh.contourNodeInds = { std::move(contourNodeInds) };

    return mesh;
}

ResultMesh construct1dMesh(const FileInput& fileInput) {
    switch (fileInput.type) {
    case 1:
        return linElem1dMesh(fileInput);
        break;

    case 2:
        return quadElem1dMesh(fileInput);
        break;

    default:
        std::cerr << "Not expected type is given to 1d mesh construct function!";
        return ResultMesh();
        break;
    }
}

ResultMesh construct2dMesh(const FileInput& fileInput) {
    switch (fileInput.type) {
    case 1:
        return tetraElem2dMesh(fileInput);
        break;

    case 2:
        return diagElem2dMesh(fileInput, false);
        break;

    case 3:
        return diagElem2dMesh(fileInput, true);
        break;

    case 4:
        return crossElem2dMesh(fileInput);
        break;

    default:
        std::cerr << "Not expected type is given to 2d mesh construct function!";
        return ResultMesh();
        break;
    }
}

void writeFile(std::ostream& file, const ResultMesh& mesh) {
    // В выводе ко всем индексам прибавляются 1, т.к. в примерах все элементы нумеруются начиная с 1

    file << mesh.NE << ' ' << mesh.NP << ' ' << mesh.NC << '\n';

    // Вывод элементов
    for (int i = 0; i < mesh.NE; ++i) {
        file << i + 1 << ' ' << mesh.elemNodeInds[i].size();
        for (int EP : mesh.elemNodeInds[i]) {
            file << ' ' << EP + 1;
        }
        file << '\n';
    }

    // Вывод узлов
    for (int i = 0; i < mesh.NP; ++i) {
        file << i + 1 << ' ' << mesh.nodes[i].x << ' ' << mesh.nodes[i].y << ' ' << mesh.nodes[i].z << '\n';
    }

    // Вывод контуров
    for (int i = 0; i < mesh.NC; ++i) {
        file << mesh.contourNodeInds[i].size() << ' ';
    }

    for (auto& contour : mesh.contourNodeInds) {
        for (int CP : contour) {
            file << '\n' << CP + 1;
        }
    }

    file << std::endl;
}

int main() {
    //std::ifstream file("C:\\Users\\Alex\\Documents\\Visual Studio Projects\\BMSTU\\Semester 7\\Programm Complex Designing\\PCD_Lab_3\\PCD_Lab_3\\test.txt");
    /*
    while (file) {
        std::string s;
        file >> s;
        std::cout << s << ' ';
    }
    */
    std::ifstream file;
    /*
#ifdef USE_INPUT_1
    file.open(input1FileName);
#else
    file.open(input2FileName);
#endif // USE_INPUT_1
    */

    file.open("input" + std::to_string(INPUT_FILE_ID) + ".txt");

    auto fileInputPtr = readFile(file);

    ResultMesh mesh;
    if (fileInputPtr->NP == 2) {
        mesh = construct1dMesh(*fileInputPtr);
    } else if (fileInputPtr->NP == 4) {
        mesh = construct2dMesh(*fileInputPtr);
    } else {
        std::cerr << "Not expected NP is given!\n";
    }

    std::ofstream outputFile(outputFileName);

    writeFile(outputFile, mesh);

    return 0;
}