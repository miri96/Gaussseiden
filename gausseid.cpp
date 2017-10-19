//Gauss - Seidel method program V 1.0.0 / Programa del metodo Gauss - Seidel V 1.0.0
#include "gausseid.h"
#include "ui_gausseid.h"

GausSeid::GausSeid(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::GausSeid)
{
    ui->setupUi(this);
    setWindowTitle("Metodo de Gauss - Seidel");// Prints the method's name at window title.
    //Connections between push buttons and slots:
    connect(ui->Inicia,SIGNAL(clicked()),this,SLOT(ejecutar()));// Button and slot for start the method.
    connect(ui->Limpia,SIGNAL(clicked()),this,SLOT(limpiar()));// Button and slot for clear the fields from user interface.
    connect(ui->Agrega,SIGNAL(clicked(bool)),this,SLOT(agregar()));// Button for add items to the table.
    // Connections between spinboxes value changed signal and the "agregar" slot:
    connect(ui->errora,SIGNAL(valueChanged(double)),this,SLOT(agregar()));
    connect(ui->imax,SIGNAL(valueChanged(int)),this,SLOT(agregar()));
    ui->Inicia->setDisabled(true);// Disables the start button till all items from table been added.
    ui->Limpia->setDisabled(true);// Disables the clear button.
    ui->errora->setDisabled(true);// Disable the accepted error spinbox till all items from table been added.
    ui->imax->setDisabled(true);// Disable the iterations spinbox till all items from table been added.
    i=0, j=0;// The row and colum counters starts from zero.
}
void GausSeid::agregar(){// With this slot the user will add items from spinboxes to the table.
    if(i<4){// From i=0 to i=3, the user can add data rows.
        ui->Limpia->setDisabled(false);// The clear button activates when the user adds a row of data.
        // All the spinboxes collects the equations variables and saves it in a 4x4 matrix.
        matrix[i][0]=ui->valorx1->value();
        matrix[i][1]=ui->valorx2->value();
        matrix[i][2]=ui->valorx3->value();
        matrix[i][3]=ui->valorx4->value();
        // This spinbox collects the value of the independent variables and saves it in a 4x1 matrix.
        indep[i][0]=ui->valori->value();
        // The introduced values shows at his respective table position.
        ui->Tabla->setItem(i,0,new QTableWidgetItem(ui->valorx1->text()));// x1 column.
        ui->Tabla->setItem(i,1,new QTableWidgetItem(ui->valorx2->text()));// x2 column.
        ui->Tabla->setItem(i,2,new QTableWidgetItem(ui->valorx3->text()));// x3 column.
        ui->Tabla->setItem(i,3,new QTableWidgetItem(ui->valorx4->text()));// x4 column.
        ui->Tabla->setItem(i,4,new QTableWidgetItem(ui->valori->text()));// = column.
        i++;// Sums one number to row counter.
    }
    if(i==4){// When i==4, the user can't add more data.
        ui->Agrega->setDisabled(true);// Disables the add button.
        ui->errora->setDisabled(false);// Enable the accepted error spinbox.
        ui->imax->setDisabled(false);// Enable the iterations spinbox.
        errorac=ui->errora->value();// Collects the value of accepted error provided by the user and saves it on "errorac".
        itera=ui->imax->value();// Collects the iterations maximum provided by user, and saves it on "itera".
    }
    if((errorac==0)||(itera==0)){// If one of the spinboxes gots a zero:
        ui->Inicia->setDisabled(true);// Disables the start button.
    }
    else{// If both spinboxes are different of zero:
        ui->Inicia->setDisabled(false);// Enables the start button.
    }
}
void GausSeid::ejecutar(){// This slot starts the method, when the start button been clicked.
    for(i=0;i<4;i++){// This cycle analyzes if the main diagonal of the main matrix have zeros.
        if(matrix[i][i]==0){// If it haves zeros:
            for(j=0;j<4;j++){// This cycle exchanges row positions for the x variables.
                if(i==3){// If is the last row.
                    ind1=indep[i][0];// A temporal variable takes the value of last row independent variable.
                    ind2=indep[i-1][0];// Another temporal variable takes the value of the previous row independent variable.
                    indep[i-1][0]=ind1;// Puts the secound independent variable on next row.
                    indep[i][0]=ind2;// Puts the first independent variable on previous row.
                    tempo[0][j]=matrix[i][j];// A temporal matrix saves the last row.
                    tempo2[0][j]=matrix[i-1][j];// Another temporal matrix saves the previous row.
                    matrix[i-1][j]=tempo[0][j];// Puts the last row on the previous row.
                    matrix[i][j]=tempo2[0][j];// Puts the previos row on the last row.
                }
                else{// If isn't the last row:
                    ind1=indep[i][0];// A temporal variable takes the value of the detected row independent variable.
                    ind2=indep[i+1][0];// Another temporal variable takes the value of the next independent variable.
                    indep[i+1][0]=ind1;// Puts the detected row independent variable on next row.
                    indep[i][0]=ind2;// Puts the next variable on previous row.
                    tempo[0][j]=matrix[i][j];// A temporal matrix saves the detected row.
                    tempo2[0][j]=matrix[i+1][j];// Another temporal matrix saves the next row.
                    matrix[i+1][j]=tempo[0][j];// Puts the detected row on the next row.
                    matrix[i][j]=tempo2[0][j];// Puts the next row on the previous row.
                }
            }
        }
    }
    double x=0, y=0, z=0, w=0;// Creates temporal variables, for save previous calculated values on method.
    double ex1=100, ex2=100, ex3=100, ex4=100;// Creates temporal variables ,for independent variables error.
    // The method starts with a cycle, and will stop when the independent calculated errors been small than the accepted error provided by the user.
    // Or when the iterations maximum been reached.
    while((ex1>errorac)&&(ex2>errorac)&&(ex3>errorac)&&(ex4>errorac)||(icont<itera)){
        x1=(indep[0][0]-(matrix[0][1]*x2)-(matrix[0][2]*x3)-(matrix[0][3]*x4))/(matrix[0][0]);// Calulates the value for x1.
        ex1=100*std::abs(x1-x)/(x1);// Calculates the error for the found value.
        x=x1;// A temporal variable saves the found value, it will be use in the next error calculation.
        x2=(indep[1][0]-(matrix[1][0]*x1)-(matrix[1][2]*x3)-(matrix[1][3]*x4))/(matrix[1][1]);// Calulates the value for x2.
        ex2=100*std::abs(x2-y)/(x2);// Calculates the error for the found value.
        y=x2;// A temporal variable saves the found value, it will be use in the next error calculation.
        x3=(indep[2][0]-(matrix[2][0]*x1)-(matrix[2][1]*x2)-(matrix[2][3]*x4))/(matrix[2][2]);// Calulates the value for x3.
        ex3=100*std::abs(x3-z)/(x3);// Calculates the error for the found value.
        z=x3;// A temporal variable saves the found value, it will be use in the next error calculation.
        x4=(indep[3][0]-(matrix[3][0]*x1)-(matrix[3][1]*x2)-(matrix[3][2]*x3))/(matrix[3][3]);// Calulates the value for x4.
        ex4=100*std::abs(x4-w)/(x4);// Calculates the error for the found value.
        w=x4;// A temporal variable saves the found value, it will be use in the next error calculation.
        icont++;// Sums one number to the iterations counter.
    } // The method finishes.
    ui->show1->display(x);// Show the last calculated value for x1;
    ui->show2->display(y);// Show the last calculated value for x2;
    ui->show3->display(z);// Show the last calculated value for x3;
    ui->show4->display(w);// Show the last calculated value for x4;
    ui->showi->display(icont);// Show the used iterations number.
}
void GausSeid::limpiar(){// This slot clears all the user interface fields, when the clear button been clicked.
    //Clears all spinboxes.
    ui->valorx1->setValue(0.0);
    ui->valorx2->setValue(0.0);
    ui->valorx3->setValue(0.0);
    ui->valorx4->setValue(0.0);
    ui->valori->setValue(0.0);
    ui->errora->setValue(0.0);
    ui->imax->setValue(0.0);
    ui->Tabla->clearContents();// Clears the table widget contents.
    // Clears all the LCD numbers.
    ui->show1->display(0.0);
    ui->show2->display(0.0);
    ui->show3->display(0.0);
    ui->show4->display(0.0);
    ui->showi->display(0.0);
    ui->Inicia->setDisabled(true);// Disables the start button.
    ui->Agrega->setDisabled(false);// Enables the add button.
    ui->imax->setDisabled(true);// The row counter goes to zero, for start over the data introduction.
    i=0;
}

GausSeid::~GausSeid()
{
    delete ui;
}
