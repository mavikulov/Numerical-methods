import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from tasks import *
from tkinter import *
from tkinter import ttk

matplotlib.use('TkAgg')


class GUI:
    def __init__(self, master):
        self.master = master
        master.title('Сплайны, команда №5, группа 382003_1')
        master.geometry('1900x1000')

        self.lf_functions = LabelFrame(master, text='Функция', width=250, height=250)
        self.lf_conditions = LabelFrame(master, text='Граничные условия', width=250, height=80)
        self.lf_functions.place(x=0, y=0)
        self.lf_conditions.place(x=307, y=0)
        self.label_conditions = Label(self.lf_conditions, text=f"S''(a) = \nS''(b) = ", background='white',
                                      font=('Consolas', 15))
        self.value_mu1 = Entry(self.lf_conditions, width=15)
        self.value_mu1.place(x=110, y=5)
        self.value_mu2 = Entry(self.lf_conditions, width=15)
        self.value_mu2.place(x=110, y=30)
        self.label_conditions.place(x=0, y=0)
        self.list_box = Listbox(self.lf_functions, selectmode=SINGLE, height=3, width=50)
        self.funcs = ['F(x) = φ(x)', 'F(x) = (ln(x+1)/x), 2 <= x <= 4',
                      'F(x) = (ln(x+1)/x + cos(10x)), 2 <= x <= 4']
        for i, func in enumerate(self.funcs):
            self.list_box.insert(i, func)
        self.list_box.pack()
        self.graphics = LabelFrame(master, text='Строить графики', width=150, height=100)
        self.graphics.place(x=0, y=70)
        self.var = IntVar()
        self.var.set(0)
        self.func_btn = Radiobutton(self.graphics, text='Функций', value=0, variable=self.var)
        self.derivative_btn = Radiobutton(self.graphics, text='Первых производных', value=1, variable=self.var)
        self.second_derivative_btn = Radiobutton(self.graphics, text='Вторых производных', value=2, variable=self.var)

        self.func_btn.pack()
        self.derivative_btn.pack()
        self.second_derivative_btn.pack()

        self.label_splits = Label(master, text='Число разбиений:')
        self.label_splits.place(x=160, y=90)
        self.entry_splits = Entry(master, background='#DCDCDC')
        self.entry_splits.place(x=270, y=90)
        self.button = Button(master, text='Интерполировать', bd=2, width=20, height=2, command=self.plot_functions)
        self.button.place(x=200, y=125)

        self.heads_coefficients = ['i', 'Xi-1', 'Xi', 'ai', 'bi', 'ci', 'di']

        self.heads_functions = ['j', "Xj", 'F(xj)', 'S(xj)', 'F(xj)-S(xj)', "S'(xj)", "F'(xj)", "F'(xj)-S'(xj)"]
        self.table_functions = ttk.Treeview(self.master, show='headings', columns=self.heads_functions, height=12)
        for header in self.heads_functions:
            self.table_functions.heading(header, text=header, anchor='center')
            if header == 'j':
                self.table_functions.column(header, anchor='center', width=50)
            else:
                self.table_functions.column(header, anchor='center', width=107)

        self.fig, self.ax = plt.subplots(1, 1)
        self.fig.set_size_inches(11, 8)
        self.ax.grid(True)
        chart = FigureCanvasTkAgg(self.fig, master)
        NavigationToolbar2Tk(chart).place(x=1530, y=800)
        chart.get_tk_widget().place(x=850, y=0)

    def clear_tables(self):
        for item in self.table_1.get_children():
            self.table_1.delete(item)
        for item in self.table_functions_1.get_children():
            self.table_functions_1.delete(item)

    def plot_functions(self):
        func = self.list_box.get(self.list_box.curselection())
        tasks = [TestTask(int(self.entry_splits.get()), -1, 1, float(self.value_mu1.get()), float(self.value_mu2.get())),
                 FirstTask(int(self.entry_splits.get()), 2, 4, float(self.value_mu1.get()), float(self.value_mu2.get())),
                 SecondTask(int(self.entry_splits.get()), 2, 4, float(self.value_mu1.get()), float(self.value_mu2.get()))]
        option = self.var.get()

        if func == 'F(x) = φ(x)':
            self.ax.clear()
            self.ax.grid(True)
            task = tasks[0]
            task()

            if option == 0:
                self.ax.plot(task.control_x, task.control_f, color='red', linestyle='--', label='F(x)')
                self.ax.plot(task.control_x, task.spline, color='blue', linestyle='--', label='S(x)')
                self.ax.legend()

                self.table_1 = ttk.Treeview(self.master, show='headings', columns=self.heads_coefficients, height=12)
                for header in self.heads_coefficients:
                    self.table_1.heading(header, text=header, anchor='center')
                    if header == 'i':
                        self.table_1.column(header, anchor='center', width=50)
                    else:
                        self.table_1.column(header, anchor='center', width=107)
                self.scrollbary_1 = Scrollbar(self.master, orient=VERTICAL)
                self.table_1.configure(yscrollcommand=self.scrollbary_1.set)
                self.scrollbary_1.configure(command=self.table_1.yview)
                self.scrollbary_1.place(x=0, y=200, width=15, height=275)
                for i in range(1, task.n + 1):
                    self.table_1.insert("", END, values=(
                        i, task.x[i - 1], task.x[i], task.a[i], task.b[i], task.c[i], task.d[i]))
                self.table_1.place(x=0, y=200)

                self.table_functions_1 = ttk.Treeview(self.master, show='headings', columns=self.heads_functions,
                                                      height=12)
                for header in self.heads_functions:
                    self.table_functions_1.heading(header, text=header, anchor='center')
                    if header == 'j':
                        self.table_functions_1.column(header, anchor='center', width=50)
                    else:
                        self.table_functions_1.column(header, anchor='center', width=113)

                ### ВТОРАЯ ТАБЛИЦА ###
                for i in range(len(task.control_x)):
                    self.table_functions_1.insert("", END, values=(i, task.control_x[i], task.control_f[i],
                                                                   task.spline[i], task.control_f[i] - task.spline[i],
                                                                   task.control_f_derivative[i],task.spline_derivative[i],
                                                                   task.control_f_derivative[i] - task.spline_derivative[i]))

                self.table_functions_1.place(x=0, y=600)
                self.scrollbary_f1 = Scrollbar(self.master, orient=VERTICAL)
                self.table_functions_1.configure(yscrollcommand=self.scrollbary_f1.set)
                self.scrollbary_f1.configure(command=self.table_functions_1.yview)
                self.scrollbary_f1.place(x=0, y=600, width=15, height=275)

                self.label_error = Label(self.master,
                                         text=f"max|f(x) - S(x)| = {np.abs(np.max(task.control_f - task.spline))} при x = {task.control_x[np.argmax(((task.control_f - task.spline)))]}\n"
                                              f"max|f'(x) - S'(x)| = {np.abs(np.max(task.control_f_derivative - task.spline_derivative))} "
                                              f"при x = {task.control_x[np.argmax(((task.control_f_derivative - task.spline_derivative)))]}\n"
                                              f"max|f''(x)- S''(x)| = {np.abs(np.max(task.control_f_second_derivative - task.spline_second_derivative))} "
                                              f"при x = {task.control_x[np.argmax(((task.control_f_second_derivative - task.spline_second_derivative)))]}",
                                         font=('Consolas', 14))
                self.label_error.place(x=0, y=500)

            elif option == 1:
                self.ax.plot(task.control_x, task.control_f_derivative, color='green', linestyle='--', label="F'(x)")
                self.ax.plot(task.control_x, task.spline_derivative, color='blue', linestyle='--', label="S'(x)")
                self.ax.legend()

            else:
                self.ax.plot(task.control_x, task.control_f_second_derivative, color='green', linestyle='--', label="F''(x)")
                self.ax.plot(task.control_x, task.spline_second_derivative, color='blue', linestyle='--', label="S''(x)")
                self.ax.legend()

            chart = FigureCanvasTkAgg(self.fig, self.master)
            NavigationToolbar2Tk(chart).place(x=1530, y=800)
            chart.get_tk_widget().place(x=850, y=0)

        elif func == 'F(x) = (ln(x+1)/x), 2 <= x <= 4':
            self.ax.clear()
            self.ax.grid(True)
            task = tasks[1]
            task()

            if option == 0:
                self.ax.plot(task.control_x, task.control_f, color='red', linestyle='--', label='F(x)')
                self.ax.plot(task.control_x, task.spline, color='blue', linestyle='--', label='S(x)')
                self.ax.legend()

                self.clear_tables()
                for i in range(1, task.n + 1):
                    self.table_1.insert("", END, values=(
                        i, task.x[i - 1], task.x[i], task.a[i], task.b[i], task.c[i], task.d[i]))
                self.table_1.place(x=0, y=200)

                for i in range(len(task.control_x)):
                    self.table_functions_1.insert("", END, values=(i, task.control_x[i], task.control_f[i],
                                                                   task.spline[i], task.control_f[i] - task.spline[i],
                                                                   task.control_f_derivative[i],task.spline_derivative[i],
                                                                   task.control_f_derivative[i] - task.spline_derivative[i]))

                self.label_error = Label(self.master,
                                         text=f"max|f(x) - S(x)| = {np.abs(np.max(task.control_f - task.spline))} при x = {task.control_x[np.argmax(((task.control_f - task.spline)))]}\n"
                                              f"max|f'(x) - S'(x)| = {np.abs(np.max(task.control_f_derivative - task.spline_derivative))} "
                                              f"при x = {task.control_x[np.argmax(((task.control_f_derivative - task.spline_derivative)))]}\n"
                                              f"max|f''(x)- S''(x)| = {np.abs(np.max(task.control_f_second_derivative - task.spline_second_derivative))} "
                                              f"при x = {task.control_x[np.argmax(((task.control_f_second_derivative - task.spline_second_derivative)))]}",
                                         font=('Consolas', 14))
                self.label_error.place(x=0, y=500)

            elif option == 1:
                self.ax.plot(task.control_x, task.control_f_derivative, color='green', linestyle='--', label="F'(x)")
                self.ax.plot(task.control_x, task.spline_derivative, color='blue', linestyle='--', label="S'(x)")
                self.ax.legend()

            else:
                self.ax.plot(task.control_x, task.control_f_second_derivative, color='green', linestyle='--', label="F''(x)")
                self.ax.plot(task.control_x, task.spline_second_derivative, color='blue', linestyle='--', label="S''(x)")
                self.ax.legend()

            chart = FigureCanvasTkAgg(self.fig, self.master)
            NavigationToolbar2Tk(chart).place(x=1530, y=800)
            chart.get_tk_widget().place(x=850, y=0)

        else:
            self.ax.clear()
            self.ax.grid(True)
            task = tasks[2]
            task()

            if option == 0:
                self.ax.plot(task.control_x, task.control_f, color='red', linestyle='--', label='F(x)')
                self.ax.plot(task.control_x, task.spline, color='blue', linestyle='--', label='S(x)')
                self.ax.legend()

                self.clear_tables()
                for i in range(1, task.n + 1):
                    self.table_1.insert("", END, values=(
                        i, task.x[i - 1], task.x[i], task.a[i], task.b[i], task.c[i], task.d[i]))
                self.table_1.place(x=0, y=200)

                for i in range(len(task.control_x)):
                    self.table_functions_1.insert("", END, values=(i, task.control_x[i], task.control_f[i],
                                                                   task.spline[i], task.control_f[i] - task.spline[i],
                                                                   task.control_f_derivative[i],
                                                                   task.spline_derivative[i],
                                                                   task.control_f_derivative[i] -
                                                                   task.spline_derivative[i]))

                self.label_error = Label(self.master,
                                         text=f"max|f(x) - S(x)| = {np.abs(np.max(task.control_f - task.spline))} при x = {task.control_x[np.argmax(((task.control_f - task.spline)))]}\n"
                                              f"max|f'(x) - S'(x)| = {np.abs(np.max(task.control_f_derivative - task.spline_derivative))} "
                                              f"при x = {task.control_x[np.argmax(((task.control_f_derivative - task.spline_derivative)))]}\n"
                                              f"max|f''(x)- S''(x)| = {np.abs(np.max(task.control_f_second_derivative - task.spline_second_derivative))} "
                                              f"при x = {task.control_x[np.argmax(((task.control_f_second_derivative - task.spline_second_derivative)))]}",
                                         font=('Consolas', 14))
                self.label_error.place(x=0, y=500)

            elif option == 1:
                self.ax.plot(task.control_x, task.control_f_derivative, color='green', linestyle='--', label="F'(x)")
                self.ax.plot(task.control_x, task.spline_derivative, color='blue', linestyle='--', label="S'(x)")
                self.ax.legend()

            else:
                self.ax.plot(task.control_x, task.control_f_second_derivative, color='green', linestyle='--', label="F''(x)")
                self.ax.plot(task.control_x, task.spline_second_derivative, color='blue', linestyle='--', label="S''(x)")
                self.ax.legend()

            chart = FigureCanvasTkAgg(self.fig, self.master)
            NavigationToolbar2Tk(chart).place(x=1530, y=800)
            chart.get_tk_widget().place(x=850, y=0)


root = Tk()
gui = GUI(root)
root.mainloop()
