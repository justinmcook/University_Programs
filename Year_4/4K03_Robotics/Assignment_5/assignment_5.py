from math import radians, sin, cos
from numpy.linalg import inv
from numpy import array, transpose

# gravity
g = -9.81

def latex_print1d(matrix):
    for i in range(matrix.shape[0]-1):
        print(str(matrix[i]) + " & ",end='')
    print(str(matrix[matrix.shape[0]-1]) + "\\\\")

def latex_print2d(matrix):
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]-1):
            print(str(matrix[i][j]) + " & ", end="")
        print(str(matrix[i][matrix.shape[1]-1]) + "\\\\")

def calculate_jacobi(theta1,theta2,theta3,theta12,theta123,a1,a2,a3):
    return array((
        [-a1*sin(theta1)-a2*sin(theta12)-a3*sin(theta123), -a2*sin(theta12)-a3*sin(theta123), -a3*sin(theta123)],
        [a1*cos(theta1)+a2*cos(theta12)+a3*cos(theta123), a2*cos(theta12)+a3*cos(theta123), a3*cos(theta123)],
        [1,1,1]
    ))

def calculate_gravity(theta1,theta2,theta3,theta12,theta123,a1,a2,a3,m1,m2,m3):
    return array((
        [(0.5*m1+m2+m3)*g*a1*cos(theta1)+(0.5*m2+m3)*g*a2*cos(theta12)+0.5*m3*g*a3*cos(theta123)],
        [(0.5*m2+m3)*g*a2*cos(theta12)+0.5*m3*g*a3*cos(theta123)],
        [0.5*m3*g*a3*cos(theta123)]
    ))

def calculate_torques(theta1,theta2,theta3,a1,a2,a3,m1,m2,m3,fWorld):
    theta12 = theta1 + theta2
    theta123 = theta1 + theta2 + theta3

    jacobi = calculate_jacobi(theta1,theta2,theta3,theta12,theta123,a1,a2,a3)
    jacobiTrans = transpose(jacobi)
    gravity = calculate_gravity(theta1,theta2,theta3,theta12,theta123,a1,a2,a3,m1,m2,m3)
    print("Part A Jacobian is:")
    latex_print2d(jacobi)
    print("Part A Jacobian transpose is:")
    latex_print2d(jacobiTrans)
    print("Part A Gravity terms is:")
    latex_print2d(gravity)
    return jacobiTrans.dot(fWorld) + gravity

def calculate_force_output(theta1,theta2,theta3,a1,a2,a3,m1,m2,m3,torques):
    theta12 = theta1 + theta2
    theta123 = theta1 + theta2 + theta3

    jacobi = calculate_jacobi(theta1,theta2,theta3,theta12,theta123,a1,a2,a3)
    jacobiTrans = transpose(jacobi)
    jacobiTransInv = inv(jacobiTrans)
    gravity = calculate_gravity(theta1,theta2,theta3,theta12,theta123,a1,a2,a3,m1,m2,m3)
    print("Part B Jacobian is:")
    latex_print2d(jacobi)
    print("Part B Jacobian transpose is:")
    latex_print2d(jacobiTrans)
    print("Part B Jacobian transpose inverse is:")
    latex_print2d(jacobiTransInv)
    print("Part B Gravity terms is:")
    latex_print2d(gravity)
    return jacobiTransInv.dot(torques-gravity)

def main():
    a1 = 0.5
    a2 = 0.5
    a3 = 0.1
    m1 = 10.0
    m2 = 10.0
    m3 = 2.0
    mPayload = 5.0

    ##### PART A #####
    theta1A = radians(45.0)
    theta2A = radians(-75.0)
    theta3A = radians(30.0)
    fWorldA = array([[0],[mPayload*g],[0]])
    # Find the torques
    torquesA = calculate_torques(theta1A,theta2A,theta3A,a1,a2,a3,m1,m2,m3,fWorldA)
    # Print the torques
    print("Part A torques are:")
    latex_print2d(torquesA)

    ##### PART B #####
    resolution = 0.1
    torques = array(([resolution],[resolution],[resolution]))
    # First configuration is same thetas as part A
    theta1B = radians(45.0)
    theta2B = radians(-5.0)
    theta3B = radians(-40.0)
    forcesA = calculate_force_output(theta1A,theta2A,theta3A,a1,a2,a3,m1,m2,m3,torques)
    forcesB = calculate_force_output(theta1B,theta2B,theta3B,a1,a2,a3,m1,m2,m3,torques)
    print("Part B configuration A forces are:")
    latex_print2d(forcesA)
    print("Part B configuration B forces are:")
    latex_print2d(forcesB)

main()