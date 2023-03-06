from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tkinter import *
from tkinter import ttk
import matplotlib
import matplotlib.pyplot as plt
import PIL
import math
import RK34


matplotlib.use('TkAgg')
window = Tk()
window.geometry('1800x950')
window.title("Численный анализ задания №9, Викулов М. команда №5")
Label(window, text='a1: ', font=('Consolas', 10)).place(x=0, y=420)
Label(window, text='a2: ', font=('Consolas', 10)).place(x=0, y=440)
Label(window, text='m: ', font=('Consolas', 10)).place(x=0, y=460)
Label(window, text='h0: ', font=('Consolas', 10)).place(x=230, y=420)
Label(window, text='u0: ', font=('Consolas', 10)).place(x=230, y=440)
Label(window, text='eps: ', font=('Consolas', 10)).place(x=230, y=460)
Label(window, text='Nmax: ', font=('Consolas', 10)).place(x=230, y=480)
Label(window, text='stop: ', font=('Consolas', 10)).place(x=230, y=500)
Label(window, text='εгр: ', font=('Consolas', 10)).place(x=230, y=520)
Label(window, text='Параметры системы', font=('Consolas', 11)).place(x=0, y=396)
Label(window, text='Условия задачи', font=('Consolas', 11)).place(x=230, y=396)
Label(window, text='Решить численно задачу Коши для ОДУ 1-го порядка:', font=('Consolas', 11)).place(x=0, y=0)
Label(window, text='Точное решение данной задачи Коши', font=('Consolas', 14)).place(x=1340, y=20)

value_a1 = Entry(window, background='#DCDCDC', width=15)
value_a2 = Entry(window, background='#DCDCDC', width=15)
value_m = Entry(window, background='#DCDCDC', width=15)
value_u0 = Entry(window, background='#DCDCDC', width=15)
value_h0 = Entry(window, background='#DCDCDC', width=15)
value_eps = Entry(window, background='#DCDCDC', width=15)
value_n_max = Entry(window, background='#DCDCDC', width=15)
value_stop = Entry(window, background='#DCDCDC', width=15)
value_bound = Entry(window, background='#DCDCDC', width=15)

value_a1.place(x=30, y=420)
value_a2.place(x=30, y=440)
value_m.place(x=30, y=460)
value_h0.place(x=270, y=420)
value_u0.place(x=270, y=440)
value_eps.place(x=270, y=460)
value_n_max.place(x=270, y=480)
value_stop.place(x=270, y=500)
value_bound.place(x=270, y=520)

# Вставка уравнения и его точного решения в виде математических выражений
canvas_equation = Canvas(window, width=400, height=100)
canvas_equation.place(x=0, y=20)
equation = PIL.ImageTk.PhotoImage(PIL.Image.open("equation№1.png"))
canvas_equation.create_image(0, 10, anchor=NW, image=equation)
canvas_solution = Canvas(window, width=1000, height=200)
canvas_solution.place(x=1300, y=50)
solution = PIL.ImageTk.PhotoImage(PIL.Image.open("solution.jpg"))
canvas_solution.create_image(0, 10, anchor=NW, image=solution)

# Проверка на нажатие галочки
checked_var = BooleanVar()
show_legend = True
show_reference = True


# Вызов функции кнопкой
def get_solution():
    global a1, a2, m
    a1 = float(value_a1.get())
    a2 = float(value_a2.get())
    m = float(value_m.get())
    h0 = float(value_h0.get())
    assert all([a1 > 0, a2 > 0, m > 0, h0 > 0]), 'Параметры a1, a2, m, h0 > 0'
    u0 = float(value_u0.get())
    n_max = int(value_n_max.get())
    x0 = 0
    h = [h0, ]
    x = [x0, ]
    v = [u0, ]
    stop = float(value_stop.get())
    assert stop >= x0, 'Значение stop >= x0'
    eps = float(value_eps.get())
    bound = float(value_bound.get())
    if checked_var.get():
        x, v, h, V2i, OLP, doublings, divisions, u = RK34.RK34_with_control(x, v, h, 0, eps, n_max, stop, a1, a2, m, u0, bound)
        difference = [vi - v2i for vi, v2i in zip(v, V2i)]
        doublings[-1] = divisions[-1] = 0
    else:
        x, v, h, u = RK34.RK34_without_control(x, v, h, 0, n_max, stop, a1, a2, m, u0, bound)
    h.insert(0, 0)

    x = [round(elem, 7) for elem in x]
    difference_between_true_and_numeric = [math.fabs(ui - vi) for ui, vi in zip(u, v)]

    # Построение таблицы
    heads = ['n', 'Xi', 'hi', 'Vi', 'V2i', 'Vi-V2i', 'Ui', 'Vi-Ui', 'ОЛП', 'C1', 'C2']
    heads_without_control = ['n', 'Xi', 'hi', 'Vi', 'V2i', 'Vi-V2i', 'Ui', 'Vi-Ui', 'C1', 'C2']
    if checked_var.get():
        table = ttk.Treeview(window, show='headings', columns=heads, height=12)
        for header in heads:
            table.heading(header, text=header, anchor='center')
            if header == "n":
                table.column(header, anchor='center', width=80)
            else:
                table.column(header, anchor='center', width=150)
        for i in range(len(x)):
            table.insert("", END,
                         values=(
                             i, x[i], h[i], v[i], V2i[i], difference[i], u[i],
                             difference_between_true_and_numeric[i],
                             OLP[i], doublings[i], divisions[i]))

        table.place(x=0, y=600)

    else:
        table_without_olp = ttk.Treeview(window, show='headings', columns=heads_without_control, height=12)
        for header in heads_without_control:
            table_without_olp.heading(header, text=header, anchor='center')
            if header == "n":
                table_without_olp.column(header, anchor='center', width=78)
            else:
                table_without_olp.column(header, anchor='center', width=167)

        for i in range(len(x)):
            table_without_olp.insert("", END,
                                     values=(
                                     i, x[i], h[i], v[i], 0, 0, u[i], difference_between_true_and_numeric[i], 0,
                                     0))
        table_without_olp.place(x=0, y=600)

    # Добавление полосы прокрутки
    scrollbary = Scrollbar(window, orient=VERTICAL)
    scrollbarx = Scrollbar(window, orient=HORIZONTAL)
    if checked_var.get():
        table.configure(yscrollcommand=scrollbary.set, xscrollcommand=scrollbarx.set)
        scrollbarx.configure(command=table.xview)
        scrollbary.configure(command=table.yview)
    else:
        table_without_olp.configure(yscrollcommand=scrollbary.set, xscrollcommand=scrollbarx.set)
        scrollbarx.configure(command=table_without_olp.xview)
        scrollbary.configure(command=table_without_olp.yview)
    scrollbary.place(x=0, y=600, width=15, height=275)

    # Построение графика решения
    global show_legend
    ax.plot(x, v, linestyle='--', label='Численное решение', color='red')
    ax.plot(x, u, linestyle='--', label='Истинное решение', color='green')
    ax.set_xlabel('x, время')
    ax.set_ylabel('u(x), скорость')
    if show_legend:
        ax.legend()
        show_legend = False
    chart = FigureCanvasTkAgg(fig, window)
    NavigationToolbar2Tk(chart).place(x=1230, y=450)
    chart.get_tk_widget().place(x=420, y=0)

    # Создание справки
    ref = LabelFrame(window, text='Справка', width=480, height=380)
    global show_reference
    ref.place(x=1240, y=160)
    if checked_var.get():
        info_with_control = Label(ref, text=f'Параметр a1 = {a1}\n'
                                            f'Параметр a2 = {a2}\n'
                                            f'Масса точки m = {m}\n'
                                            f'В конечный момент времени x = {x[-1]} значение скорости u(x) = {v[-1]}   \n'
                                            f'Начальный момент времени x0 = 0\n'
                                            f'Начальное значение скорости u(0) = {u0}\n\n'
                                            f'Максимальный шаг = {max(h[:len(h) - 1])} при x = {round(x[h.index(max(h[:len(h) - 1]))], 6)}\n'
                                            f'Минимальный шаг = {min(h[1:])} при x = {round(x[h.index(min(h[1:]))], 6)}\n'
                                            f'Число удвоений шага = {sum(doublings[:len(doublings) - 1])}\n'
                                            f'Число делений шага = {sum(divisions[:len(divisions) - 1])}\n'
                                            f'Максимальная оценка лок. погрешности max|ОЛП| = {max(OLP[:len(OLP) - 1])}\n'
                                            f'Максимальная глобальная погрешность |E|= {max(difference_between_true_and_numeric)}',
                                  font=('Calibri', 11),
                                  width=70, anchor=NW)
        print(info_with_control['height'])
        info_with_control.pack()
    else:
        info_without_control = Label(ref, text=f'Параметр a1 = {a1}\n'
                                               f'Параметр a2 = {a2}\n'
                                               f'Масса точки m = {m}\n'
                                               f'В конечный момент времени x = {x[-1]} значение скорости u(x) = {v[-1]}   \n'
                                               f'Начальный момент времени x0 = 0\n'
                                               f'Начальное значение скорости u(0) = {u0}\n\n'
                                               f'Максимальный шаг = {h0}\n'
                                               f'Минимальный шаг = {h[-2]}\n'
                                               f'Число удвоений шага = {0}\n'
                                               f'Число делений шага = {0}\n'
                                               f'Максимальная глобальная погрешность |E| = {max(difference_between_true_and_numeric)}',
                                     font=('Calibri', 11),
                                     width=70, height=13, anchor=NW)
        info_without_control.pack()
    show_reference = False


# Описание параметров
labelframe = LabelFrame(window, text='Метод решения и описание параметров', width=100)
labelframe.place(x=0, y=130)
left = Label(labelframe, text=f"Вариант задания №3\n"
                              f"Метод Рунге-Кутта явный 3-го порядка (p = 3)\n"
                              f"u(x) - скорость точки в момент времени x\n"
                              f"u0 - начальная скорость точки в момент времени x=0\n"
                              f"R(u) - сила сопротивления, зависящая от скорости\n"
                              f"a1, a2 - постоянные (a1, a2 > 0)\n"
                              f"m - масса точки\n"
                              f"h0 - начальный шаг\n"
                              f"eps - параметр контроля лок. погрешности\n"
                              f"Nmax - максимальное допустимое количество итераций\n"
                              f"stop - конечное значение для рассчета траектории\n"
                              f"εгр - значение для контроля выхода за правую границу\n"
                              f"ОЛП - оценка локальной погрешности S* = 8 * S",
             justify=LEFT, font=("Consolas", 11))
left.pack()
fig, ax = plt.subplots(1, 1)
fig.set_size_inches(8, 5)
ax.set_title("Графики скоростей v(x), u(x) численного и истинного решений")
ax.grid(True)
chart = FigureCanvasTkAgg(fig, window)
NavigationToolbar2Tk(chart).place(x=1230, y=450)
chart.get_tk_widget().place(x=420, y=0)

button = Button(window, text='Построить решение', command=get_solution, bd=2, width=20, height=2, )
button.place(x=0, y=490)
check_button = Checkbutton(window, text='Использовать контроль локальной погрешности', variable=checked_var)
check_button.place(x=0, y=540)
window.mainloop()
