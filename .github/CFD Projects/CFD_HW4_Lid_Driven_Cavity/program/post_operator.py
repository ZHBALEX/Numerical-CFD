import numpy as np
import re
import matplotlib.pyplot as plt

def load_2d_arrays_with_numpy(label):
    filename = f'output_N={label}.txt'
    try:
        with open(filename, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print("文件未找到，请检查文件路径是否正确")
        return None, None, None

    # 辅助函数来提取并构建二维数组
    def extract_2d_array(array_content):
        array_data = []
        # 匹配方括号内的所有内容，包括可能的多行
        matches = re.findall(r'\[(.*?)]', array_content.replace('\n', ''), re.DOTALL)
        for match in matches:
            # 替换掉可能存在的换行符，然后解析数值
            formatted_match = ' '.join(match.split())
            numbers = np.fromstring(formatted_match, sep=' ')
            array_data.append(numbers)
        return np.array(array_data)

    # 改进正则表达式以匹配整个数组块
    u_data = re.search(r'Array u:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    v_data = re.search(r'Array v:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)
    p_data = re.search(r'Array p:\s*((?:\[\s*.*?\s*\][,\s]*)+)', content, re.DOTALL)

    if not (u_data and v_data and p_data):
        print("未能正确提取u, v, 或 p的数据。请检查数组标记和文件格式。")
        return np.array([]), np.array([]), np.array([])

    u = extract_2d_array(u_data.group(1))
    v = extract_2d_array(v_data.group(1))
    p = extract_2d_array(p_data.group(1))

    return u, v, p

def smalluvp(u,v,p):
    return u[1:-1,1:-1],v[1:-1,1:-1],p[1:-1,1:-1] 

def Center(u,v,p, label):
    u,v,p = smalluvp(u,v,p)
    
    center = int(label/2)
    u_x_center = (u[:,center] + u[:,center-1])/2
    v_x_center = (v[:,center] + v[:,center-1])/2
    p_x_center = (p[:,center] + p[:,center-1])/2

    u_y_center = (u[center,:] + u[center-1,:])/2
    v_y_center = (v[center,:] + v[center-1,:])/2
    p_y_center = (p[center,:] + p[center-1,:])/2

    return u_x_center, v_x_center, p_x_center\
    ,u_y_center, v_y_center, p_y_center


# def drawing(u_x_c, u_y_c):
    
    
def get_ue_x_c():
    ye_0 = [1.00000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
    ue_x_c_0 = [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000]
    return ye_0[::-1], ue_x_c_0[::-1]

def get_ue_y_c():
    xe_0 = [1.0000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]
    xe_y_c_0 = [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000]
    return xe_0[::-1], xe_y_c_0[::-1]


# def drawing(u_x_c_16, u_x_c_64, u_x_c_128):
#     ye, ue_x_c = get_ue_x_c()
#     y = np.linspace(0,N*1/N, N)
#     plt.figure(figsize=(8, 8))
#     plt.plot(ye, ue_x_c, label='Paper Result')
#     plt.plot(y, u_x_c, label='Result')
#     plt.xlabel('X')
#     plt.ylabel('u velocity')
#     plt.legend()
#     # plt.title('u Contour')
#     # plt.savefig(f'u_contour_{label}.png')
#     plt.show()


def drawing_x_c(u_x_c_16, u_x_c_64, u_x_c_128):
    ye, ue_x_c_vertical = get_ue_x_c()
    y_16 = np.linspace(0,1, 16)
    y_64 = np.linspace(0,1, 64)
    y_128 = np.linspace(0,1, 128)
    plt.figure(figsize=(8, 8))
    plt.plot( ue_x_c_vertical, ye, '-x', label='Paper Result', color= 'red')
    plt.plot( u_x_c_16,y_16, '-x', label='Result_16',color= 'blue')
    plt.plot( u_x_c_64,y_64, '-x', label='Result_64',color= 'green')
    plt.plot( u_x_c_128,y_128, '-', label='Result_128',color= 'navy')
    plt.xlabel('Y')
    plt.ylabel('u velocity')
    plt.legend()
    plt.title('u along verticle line through geometric center')
    plt.savefig('u along verticle line through geometric center')
    # plt.show()



    

def drawing_y_c(v_y_c_16, v_y_c_64, v_y_c_128):
    xe, ue_y_c_vertical = get_ue_y_c()
    plt.figure(figsize=(12,12))
    x_16 = np.linspace(0,16*1/16, 16)
    x_64 = np.linspace(0,64*1/64, 64)
    x_128 = np.linspace(0,128*1/128, 128)

    plt.plot(xe, ue_y_c_vertical, '-x', label='Paper Result',color= 'red')
    plt.plot(x_16, v_y_c_16, '-x', label='Result_16',color= 'blue')
    plt.plot(x_64, v_y_c_64, '-x', label='Result_64',color= 'green')
    plt.plot(x_128, v_y_c_128, '-', label='Result_128',color= 'navy')
    plt.xlabel('X')
    plt.ylabel('v velocity')
    plt.legend()
    plt.title('v along horizontal line through geometric center')
    plt.savefig('v along horizontal line through geometric center')
    # plt.show()


def center_get(N):
    u, v, p = load_2d_arrays_with_numpy(N)
    return Center(u,v,p, N)








def main():
    # 使用示例
    # N = 128
    # u, v, p = load_2d_arrays_with_numpy(N)
    # print(u[-2,:])
    # u_x_c, v_x_c, p_x_c, u_y_c, v_y_c, p_y_c = Center(u,v,p, N)

    u_x_c_128, v_x_c_128, p_x_c_128, u_y_c_128, v_y_c_128, p_y_c_128 = center_get(128)
    u_x_c_64, v_x_c_64, p_x_c_64, u_y_c_64, v_y_c_64, p_y_c_64 = center_get(64)
    u_x_c_16, v_x_c_16, p_x_c_16, u_y_c_16, v_y_c_16, p_y_c_16 = center_get(16)

    # drawing_x_c(u_x_c_16, u_x_c_64, u_x_c_128)
    # drawing_y_c(v_y_c_16, v_y_c_64, v_y_c_128)



if __name__ == '__main__':
    main()

    








