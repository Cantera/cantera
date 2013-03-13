from Tkinter import *

def make_menu(name, menubar, list):
    from types import *
    button=Menubutton(menubar, text=name, padx=3,pady=1)
    button.pack(side=LEFT, anchor=W)
    menu = Menu(button,tearoff=FALSE)
    for entry in list:
        if entry == 'separator':
            menu.add_separator({})
        elif type(entry)==ListType:
            for num in entry:
                menu.entryconfig(num,state=DISABLED)
        elif type(entry[1]) != ListType:
            if len(entry) == 2 or entry[2] == 'command':
                menu.add_command(label=entry[0],
                                 command=entry[1])
            elif entry[2] == 'check':
                entry[3].set(0)
                if len(entry) >= 5: val = entry[4]
                else: val = 1
                menu.add_checkbutton(label=entry[0],
                                     command=entry[1], variable = entry[3],
                                     onvalue=val)
            else:
                submenu=make_menu(entry[0],menu, entry[1])
                menu.add_cascade(label=entry[0],
                                 menu=submenu)

    button['menu']=menu
    return button
