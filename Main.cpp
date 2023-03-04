/*
**********************************************************
By Sreyash (Srey) Raychaudhuri
Delaunay Triangulation Demo
**********************************************************
* */

#include <glad/glad.h>
#include <glfw/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <iostream>
#include <vector>

#include "DelaunayTriangulator.h"

const int SCREEN_WIDTH = 1920;
const int SCREEN_HEIGHT = 1080;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window, int& currState, int& prevState, bool& triangulateNextStep);

int main() {

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(1080, 1080, "PCG Dungeon Gen.", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);


    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }


    std::unique_ptr<DelaunayTriangulator> p_Triangulator
        = std::make_unique<DelaunayTriangulator>();

    int currState, prevState;
    bool triangulateNextStep = false;

    glColor3f(1, 1, 1);

    p_Triangulator->Init(10);


    // Update loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        processInput(window, currState, prevState, triangulateNextStep);
        glClear(GL_COLOR_BUFFER_BIT);
        glClearColor(0, 0, 0, 1.0f);

        if (triangulateNextStep) {
            p_Triangulator->Triangulate();
            triangulateNextStep = false;
        }

        p_Triangulator->Display();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }


    glfwTerminate();
    return 0;
}

// To make sure the viewport matches the new window dimensions; 
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window, int& currState, int& prevState, bool& triangulateNextStep) {
    currState = glfwGetKey(window, GLFW_KEY_P);
    if (currState == GLFW_PRESS && prevState == GLFW_RELEASE) {
        triangulateNextStep = true;
    }
    prevState = currState;
}