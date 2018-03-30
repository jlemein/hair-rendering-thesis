#include <iostream>
#include "sayhello.h"

int main() {
    std::cout << "render-frame <options>" << std::endl;\

    SayHello sayHello;
    sayHello.sayHello("render-frame");
}
