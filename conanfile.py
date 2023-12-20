from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps

class PonioConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = [
        "ninja/1.11.1",
        "doctest/2.4.10",
        "eigen/3.4.0",
        "cli11/2.3.2"
    ]

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
