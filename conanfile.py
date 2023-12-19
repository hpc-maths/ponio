from conans import ConanFile, CMake


class PonioConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = [
        "doctest/2.4.10",
        "eigen/3.4.0",
        "cli11/2.3.2"
    ]
    generators = ["CMakeDeps", "CMakeToolchain"]

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
