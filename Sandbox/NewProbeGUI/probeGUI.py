import tkinter as tk
from tkinter import ttk


def create_gui():
    root = tk.Tk()
    root.title('Probe GUI')

    font = 'Cambria'
    borderwidth = 3
    relief = tk.GROOVE
    frame1 = ttk.Frame(root, relief=relief, borderwidth=borderwidth, height=800, width=400)
    frame2 = ttk.Frame(root, relief=relief, borderwidth=borderwidth, height=400, width=400)
    frame3 = ttk.Frame(root, relief=relief, borderwidth=borderwidth, height=400, width=400)
    frame4 = ttk.Frame(root, relief=relief, borderwidth=borderwidth, height=400, width=800)
    frame1.grid(row=0, column=0, rowspan=2)
    frame2.grid(row=0, column=1)
    frame3.grid(row=0, column=2)
    frame4.grid(row=1, column=1, columnspan=2)

    # Create the readout box.
    frame4.pack_propagate(False)
    textbox = tk.Text(frame4, background='black')
    textbox.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

    # Create the Checkbutton objects to hold probe selection.
    frame1.grid_propagate(False)
    ttk.Label(frame1, text='Probes to retrieve data for:', font=(font, 20, 'underline'),
              padding=7).grid(row=0, column=0, columnspan=3)
    all_A_probes = ['A2',  'A15', 'A17', 'A18', 'A19', 'A20', 'A21', 'A22',
                    'A23', 'A24', 'A25', 'A27', 'A28', 'A31', 'A32', 'A33',
                    'A34', 'A35']
    all_B_probes = ['B2', 'B7', 'B8', 'B10']
    all_C_probes = ['C2', 'C7', 'C8', 'C10']
    all_A_checks = []; all_B_checks = []; all_C_checks = []
    for p in all_A_probes:
        all_A_checks.append(ttk.Checkbutton(frame1, variable=p, text=p))
    for p in all_B_probes:
        all_B_checks.append(ttk.Checkbutton(frame1, variable=p, text=p))
    for p in all_C_probes:
        all_C_checks.append(ttk.Checkbutton(frame1, variable=p, text=p))

    # Put in rows in frame1.
    A_row_idx = 1
    for checkbutton in all_A_checks:
        checkbutton.grid(row=A_row_idx, column=0, sticky='w', padx=20, pady=3)
        A_row_idx += 1
    B_row_idx = 1
    for checkbutton in all_B_checks:
        checkbutton.grid(row=B_row_idx, column=1, sticky='w', padx=20, pady=3)
        B_row_idx += 1
    C_row_idx = 1
    for checkbutton in all_C_checks:
        checkbutton.grid(row=C_row_idx, column=2, sticky='w', padx=20, pady=3)
        C_row_idx += 1

    # At bottom add 'all probes' option.
    all_button = ttk.Checkbutton(frame1, text='Select all probes')
    all_button.grid(row=A_row_idx, column=0, columnspan=3)

    style=ttk.Style()
    style.configure('TCheckbutton', font=(font, 18))

    # Add a progress bar and logo to frame2.
    progressbar = ttk.Progressbar(frame2, orient=tk.HORIZONTAL, length=100)
    frame2.pack_propagate(False)
    progressbar.pack(fill=tk.X, expand=True, padx=10, side=tk.TOP)
    utk_logo = tk.PhotoImage(file='utk_logo.gif').subsample(4, 4)
    logo_label = ttk.Label(frame2)
    logo_label.config(image=utk_logo, borderwidth=3, relief=tk.SOLID)
    logo_label.pack(pady=25)

    # Fill in the save data box with the widgets.
    frame3.grid_propagate(False)
    ttk.Label(frame3, text='Save data as...', font=(font, 20, 'underline'), padding=7).grid(row=0, column=0, columnspan=2)
    save_text = tk.Text(frame3, width=40, height=1)
    save_text.grid(row=1, column=0, padx=5)
    file_format = tk.StringVar()
    save_opts = ttk.Combobox(frame3, textvariable=file_format, width=5)
    save_opts.grid(row=1, column=1)
    save_opts.config(values=('.csv', 'hdf5', '.mat'))
    save_button = ttk.Button(frame3, text='Save Data')
    save_button.grid(row=2, column=0, columnspan=2, pady=10)

    root.mainloop()
