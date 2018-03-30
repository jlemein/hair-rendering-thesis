#include "sayhello.h"
#include <iostream>

void SayHello::sayHello(const char* name) {
    std::cout << "Saying hello from library to " << name << std::endl;
}
