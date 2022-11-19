#ifndef UDRAWWIDGET_H
#define UDRAWWIDGET_H

#include <QWidget>
#include <vector>
#include "mainwindow.h"
#include <QPushButton>


struct vector2d
{
    long double x,y;
    vector2d() { }
    vector2d(long double ix, long double iy) : x(ix), y(iy) { }
    static double dist(vector2d A, vector2d B)
    {
        return sqrt(pow(A.x-B.x, 2) + pow(A.y-B.y, 2));
    }
    vector2d operator+(const vector2d& oth) const
    {
      return vector2d(x + oth.x, y + oth.y);
    }

    vector2d operator-(const vector2d& oth) const
    {
      return vector2d(x - oth.x, y - oth.y);
    }

    vector2d operator*(double oth) const
    {
      return vector2d(x * oth, y * oth);
    }

    vector2d operator+=(const vector2d oth)
    {
      x += oth.x;
      y += oth.y;
      return *this;
    }
     bool operator==(vector2d oth)
     {
         return (this->x == oth.x) && (this->y == oth.y);
     }
};


struct edge
{
    vector2d A, B;
    long double xx, yy;
    int flag;
    edge (vector2d vec1, vector2d vec2) : A(vec1), B(vec2)
    {
        vector2d mod = B - A;
        xx = mod.x;
        yy = mod.y;
    }
    bool operator==(edge oth)
    {
        return (A == oth.A) && (B == oth.B);
    }
    double norm()
    {
        return vector2d::dist(A, B);
    }
    double ang(vector2d vec, int flag)
    {
        edge a = edge(A, vec), b = edge(B, vec);
        if (a.norm() == 0 || b.norm() == 0 || vec == A || vec == B) return 0;
        int w = 1, buf = (xx*a.yy - yy*a.xx);
        if (buf != 0) buf = buf/abs(buf);
        if (buf == 0 || flag == buf) w = 0;
        return w*acos((pow(a.norm(), 2)+pow(b.norm(), 2)-pow(norm(),2))/(2*a.norm()*b.norm()));
    }
};

typedef enum { NoPoint, Point } pointCatched;

class UDrawWidget : public QWidget
{
    Q_OBJECT

protected:

    int pr_mouseX;
    int pr_mouseY;
    pointCatched catched_point;
public slots:
   void Clear();

public:
    explicit UDrawWidget(QWidget *parent = nullptr);

    virtual void paintEvent(QPaintEvent *event);
    virtual void mousePressEvent(QMouseEvent* pe);   // методы обработки события мыши при нажатии клавиши мыши
    virtual void mouseMoveEvent(QMouseEvent* pe);    // методы обработки события мыши при перемещении мыши
    virtual void mouseReleaseEvent(QMouseEvent* pe); // методы обработки событий мыши при отжатии клавиши мыши
    bool buffer = false;
    bool H, L;
    int N = 0;
    vector2d k;
    void Print(Ui::MainWindow *ui);
    QPushButton*a;
    UDrawWidget*b;

private:
    std::vector <vector2d> point, result;
    //std::vector <edge> result;


signals:

};

#endif // UDRAWWIDGET_H
