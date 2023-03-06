from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
import matplotlib.pyplot as plt
import RK


window = Tk()
window.geometry('1700x870')
window.title('Метод Рунге-Кутта 4 порядка для ДУ 2 порядка')

Label(window, text='Введите значение a: ').place(x=0, y=0)
Label(window, text='Введите значение b: ').place(x=0, y=20)
Label(window, text='Введите значение u01: ').place(x=0, y=40)
Label(window, text='Введите значение u02: ').place(x=0, y=60)
Label(window, text='Введите значение h0: ').place(x=0, y=80)
Label(window, text='Введите значение eps: ').place(x=0, y=100)
Label(window, text='Введите значение Nmax: ').place(x=0, y=120)
Label(window, text='Введите значение stop: ').place(x=0, y=140)

value_a = Entry(window)
value_b = Entry(window)
value_u01 = Entry(window)
value_u02 = Entry(window)
value_h0 = Entry(window)
value_eps = Entry(window)
value_n_max = Entry(window)
value_stop = Entry(window)

value_a.place(x=140, y=0)
value_b.place(x=140, y=20)
value_u01.place(x=140, y=40)
value_u02.place(x=140, y=60)
value_h0.place(x=140, y=80)
value_eps.place(x=140, y=100)
value_n_max.place(x=140, y=120)
value_stop.place(x=140, y=140)

checked_var = BooleanVar()


def get_solution():
    a = float(value_a.get())
    b = float(value_b.get())
    assert all([a > 0, b > 0, float(value_h0.get()) > 0]), 'Значения a, b, h > 0'
    v = np.array([[float(value_u01.get()), ], [float(value_u02.get()), ]])
    h = np.array([float(value_h0.get()), ])
    x0 = 0
    x = np.array([x0, ])
    eps = float(value_eps.get())
    try:
        n_max = int(value_n_max.get())
    except TypeError:
        print('Значение Nmax принимает числовые значения')
    stop = float(value_stop.get())
    assert stop >= 0, 'Значение stop >= start'
    f = (RK.f_1, RK.f_2)
    if checked_var.get():
        x, v, v2i, h, OLP, doublings, divisions = RK.RK45(x, v, h, f, x0, a, b, eps, n_max, stop, True)
    else:
        x, v, v2i, h, OLP, doublings, divisions = RK.RK45(x, v, h, f, x0, a, b, eps, n_max, stop, False)
    h = np.insert(h, [0], 0)
    x = np.around(x, decimals=10)
    v2i = np.array(v2i).transpose()
    doublings.append(0)
    divisions.append(0)

    # Построение таблицы значений
    heads = ['Xi', 'hi', 'Vi', 'V2i', 'Vi-V2i', 'OLP', 'C1', 'C2']
    table = ttk.Treeview(window, show='headings', columns=heads, height=7)
    for header in heads:
        table.heading(header, text=header, anchor='center')
        table.column(header, anchor='center')

    if checked_var.get():
        for i in range(len(x)):
            if i == len(x) - 1:
                table.insert("", END,
                             values=(x[i], h[i + 1], v[:, i], v2i[:, i], v[:, i] - v2i[:, i], OLP[i], doublings[i], divisions[i]))
            else:
                table.insert("", END,
                             values=(x[i], h[i], v[:, i], v2i[:, i], v[:, i] - v2i[:, i], OLP[i], doublings[i], divisions[i]))
    else:
        for i in range(len(x)):
            table.insert("", END, values=(x[i], h[i], v[:, i], 0, 0, 0, 0, 0))

    cols = [f'#{col}' for col in range(1, 9)]
    for number_of_col in cols:
        if number_of_col in ('#3, #4, #5, #6, #8'):
            table.column(number_of_col, stretch=NO, width=200)
        else:
            table.column(number_of_col, stretch=NO, width=60)

    # Добавление полосы прокрутки
    table.place(x=0, y=650)
    scrollbary = Scrollbar(window, orient=VERTICAL)
    scrollbarx = Scrollbar(window, orient=HORIZONTAL)
    table.configure(yscrollcommand=scrollbary.set, xscrollcommand=scrollbarx.set)
    scrollbarx.configure(command=table.xview)
    scrollbary.configure(command=table.yview)
    scrollbary.place(x=0, y=650, width=15, height=165)


    # Вывод некоторых параметров
    info_with_control = Label(window, text=f'Максимальный шаг = {np.max(h[1: h.size])}\n'
                              f'Минимальный шаг = {np.min(h[1: h.size])}\n'
                              f'Число удвоений шага = {sum(doublings)}\n'
                              f'Число делений шага = {sum(divisions)}\n'
                              f'Максимальная ошибка = {max(OLP)}')

    info_without_control = Label(window, text=f'Максимальный шаг = {value_h0.get()}\n'
                                              f'Минимальный шаг = {value_h0.get()}\n'
                                              f'Число удвоений шага = {0}\n'
                                              f'Число делений шага = {0}\n'
                                              f'Максимальная ошибка = {0}')

    if checked_var.get():
        info_with_control.place(x=1200, y=670)
    else:
        info_without_control.place(x=1200, y=670)

    # Построение графиков
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(11, 5)
    ax[0].plot(x, v[0, :], label="v1(x)", color='r')
    ax[0].set_title('График зависимости v1(x)')
    ax[0].legend()
    ax[0].grid(True)
    ax[1].plot(x, v[1, :], label='v2(x)', color='r')
    ax[1].set_title('График зависимости v2(x)')
    ax[1].legend()
    ax[1].grid(True)
    ax[2].plot(v[0, :], v[1, :], label='v2(v1)', color='b')
    ax[2].set_title('Фазовый потрет системы')
    ax[2].legend()
    ax[2].grid(True)
    chart = FigureCanvasTkAgg(fig, window)
    chart.get_tk_widget().place(x=500, y=0)


button = Button(window, text='Построить решение', command=get_solution)
button.place(x=275, y=70)
check_button = Checkbutton(window, text='Использовать контроль локальной погрешности', variable=checked_var)
check_button.place(x=0, y=200)

window.mainloop()
