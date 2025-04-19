#include <iostream>
#include "json.hpp"

int main(int argc, char *argv[])
{
    int a = std::stoi(argv[1]);
    int b = std::stoi(argv[2]);

    int result = a + b;

    nlohmann::json response;
    response["a"] = a;
    response["b"] = b;
    response["result"] = result;

    std::cout << response.dump();
    return 0;
}
