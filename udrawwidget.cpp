#include "udrawwidget.h"
#include <QPainter>
#include <QMouseEvent>
#include <QPen>
#include <QMetaObject>
#include <iostream>
#include <mainwindow.h>
#include <ui_mainwindow.h>
#include <queue>
#include <utility>

using std::vector;
using std::pair;


vector <vector2d> hermite_splain(vector <vector2d> point)
{
    int n = point.size();
    vector <vector2d> result, m(n);
    vector <double> t;
    t.push_back(0);
    for(int i = 0; i < n - 1; ++i)
        t.push_back(sqrt(pow(point[i + 1].x - point[i].x, 2) + pow(point[i + 1].y - point[i].y, 2)) + t[i]);
    double min = 0, max = t[n - 1];
    for (int i = 1; i < n - 1; ++i)
    {
        m[i].x = (((point[i + 1].x - point[i].x) / (t[i + 1] - t[i])) + ((point[i].x - point[i - 1].x) / (t[i] - t[i - 1]))) /2;
        m[i].y = (((point[i + 1].y - point[i].y) / (t[i + 1] - t[i])) + ((point[i].y - point[i - 1].y) / (t[i] - t[i - 1]))) /2;
    }
    m[0].x = (point[1].x - point[0].x) / (t[1] - t[0]);
    m[n - 1].x = (point[n - 1].x - point[n - 2].x) / (t[n - 1] - t[n - 2]);
    m[0].y = (point[1].y - point[0].y) / (t[1] - t[0]);
    m[n - 1].y = (point[n - 1].y - point[n - 2].y) / (t[n - 1] - t[n - 2]);
    for (int i = 0; min+i*1<=max; i++)
    {
        double T =  min+i*1;
        int k = 0;
        for(int j = 0; t[j] < T; ++j)
            k = j;
        double t_x = (T - t[k]) / (t[k + 1] - t[k]);
        double h00 = 2 * pow(t_x, 3) - 3 * pow(t_x, 2) + 1;
        double h10 = pow(t_x, 3) - 2 * pow(t_x, 2) + t_x;
        double h01 = -2 * pow(t_x, 3) + 3 * pow(t_x, 2);
        double h11 = pow(t_x,3) - pow(t_x,2);
        double x = h00 * point[k].x + h10 * (t[k+1] - t[k]) * m[k].x+h01*point[k+1].x + h11 * (t[k+1] - t[k]) * m[k+1].x;
        double y = h00 * point[k].y + h10 * (t[k+1] - t[k]) * m[k].y+h01*point[k+1].y + h11 * (t[k+1] - t[k]) * m[k+1].y;
        result.push_back(vector2d(x, y));
    }
    return result;
}

vector <vector2d> lagrange_splain(vector <vector2d> point)
{
    std::sort(point.begin(), point.end(),[] (const auto &x, const auto &y) { return x.x < y.x; });
    int N = point.size();
    long double min = point[0].x;
    long double max = point[N-1].x;
    vector <vector2d> result;
    for (int i = 0; min + i*1 <= max; ++i)
    {
        double x = min + i*1;
        double y = 0;
        for (int i = 0; i < N; ++i)
        {
            double Li = 1;
            for (int j = 0; j < N; ++j)
                if (i != j)
                    Li = Li * (x - point[j].x) / (point[i].x - point[j].x);
            y = y + point[i].y * Li;
        }
        result.push_back(vector2d(x, y));
    }
    return result;
}
vector<double> TDMA(vector<double> a, vector<double> b, vector<double> c, vector<double> d)
{
    int n = b.size();
    if (n < 2)
        d[0] /= b[0];
    else
    {
        n--;
        c[0] /= b[0];
        d[0] /= b[0];
        for (int i = 1; i < n; i++)
        {
            c[i] /= b[i] - a[i]*c[i-1];
            d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
        }

        d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);
        for (int i = n; i-- > 0;)
        {
            d[i] -= c[i]*d[i+1];
        }
    }
        return d;
 }

vector <vector2d> cub_splain(vector <vector2d> point)
{
    int n = point.size();
    vector <vector2d> result;
    vector<double> t;
    t.push_back(0);
    for(int i = 0; i < n - 1; ++i)
        t.push_back(sqrt(pow(point[i + 1].x - point[i].x, 2) + pow(point[i + 1].y - point[i].y, 2)) + t[i]);
    double min = 0, max = t[n - 1];
    vector<double> h(n - 1), mdx(n-1), mdy(n-1);

    for (int i = 0; i < n - 1; ++i)
        h[i] = t[i + 1] - t[i];

    for (int i = 0; i < n - 1; ++i)
    {
        mdx[i] = (point[i + 1].x - point[i].x) / h[i];
        mdy[i] = (point[i + 1].y - point[i].y) / h[i];
    }

    vector<double> a(n - 2), b(n-2), c(n-2), dx(n - 2), dy(n - 2);

    for (int i = 0; i < n - 2; ++i)
    {
        a[i] = h[i];
        c[i] = h[i+1];
        b[i] = 2 * (h[i] + h[i + 1]);
        dx[i] = 6 * (mdx[i + 1] - mdx[i]);
        dy[i] = 6 * (mdy[i + 1] - mdy[i]);
    }

    vector<double> m_x, m_y;

    if(n > 2)
    {
        m_x = TDMA(a, b, c, dx);
        m_y = TDMA(a, b, c, dy);
    }

    m_x.insert(m_x.begin(), 0);
    m_x.push_back(0);
    m_y.insert(m_y.begin(), 0);
    m_y.push_back(0);

    vector <vector2d>  s0(n - 1), s1(n - 1), s2(n - 1), s3(n - 1);

    for (int i = 0; i < n - 1; ++i)
    {
        s3[i].x = (m_x[i + 1] - m_x[i]) / (6 * h[i]);
        s2[i].x = m_x[i] / 2;
        s1[i].x = mdx[i] - h[i] * (2 * m_x[i] + m_x[i + 1]) /6;
        s0[i].x = point[i].x;

        s3[i].y = (m_y[i + 1] - m_y[i]) / (6 * h[i]);
        s2[i].y = m_y[i] / 2;
        s1[i].y = mdy[i] - h[i] * (2 * m_y[i] + m_y[i + 1]) / 6;
        s0[i].y = point[i].y;
    }

    for (int i = 0; min + i * 1 <= max; ++i)
    {
        double T = min + i * 1;
        int k = 0;
        for(int j = 0; t[j] < T; j++)
            k = j;
        double s = T - t[k];
        double x = s3[k].x * pow(s, 3) + s2[k].x * pow(s, 2) + s1[k].x * s + s0[k].x;
        double y = s3[k].y * pow(s, 3) + s2[k].y * pow(s, 2) + s1[k].y * s + s0[k].y;
        result.push_back(vector2d(x, y));
    }
    return result;
}


UDrawWidget::UDrawWidget(QWidget *parent)
    : QWidget{parent}
{
    catched_point = NoPoint;
}


void UDrawWidget::paintEvent(QPaintEvent *event)
{
    QPainter p(this);
    for (auto x : this->point)
    {
        p.setPen(Qt::red);
        p.drawEllipse(int(x.x) - 2.5, int(x.y) - 2.5, 5, 5);
    }
    p.setPen(Qt::black);
    if (point.size() > 1) result = cub_splain(point);
    if (this->result.size() > 1)
        {
            for (auto k = this->result.begin(); k != this->result.end()-1; k++)
                p.drawLine(int(k->x), int(k->y), int((k+1)->x), int((k+1)->y));
        }

}


void UDrawWidget::mousePressEvent(QMouseEvent* pe)  // методы обработки события мыши при нажатии клавиши мыши
{
    if(pe->button() == Qt::LeftButton)
    {
        bool flag = false;

        if (int(point.size()) == 0)
        {
            this->point.push_back(vector2d(pe->position().x(), pe->position().y()));
        }
        for (int i = 0; i < int(point.size()); i++)
        {
            if (vector2d::dist(point[i], vector2d(pe->position().x(), pe->position().y())) <= 2.5)
            {
                catched_point = Point;
                k=point[i];
                flag = true;
                pr_mouseX = pe->x();
                pr_mouseY = pe->y();
            }
        }
        if (!flag)
        {
            this->point.push_back(vector2d(pe->position().x(), pe->position().y()));
        }
    }
    if(pe->button() == Qt::RightButton)
    {
        auto buf = remove_if(this->point.begin(), this->point.end(), [&](vector2d a)
        {return vector2d::dist(a, vector2d(pe->position().x(), pe->position().y())) <= 2.5;});
        this->point.erase(buf, this->point.end());
    }
     update();

}

void UDrawWidget::Clear()  // очищение экрана
{
     this->point.erase(this->point.begin(), this->point.end());
     this->result.erase(this->result.begin(), this->result.end());
     update();

}

void UDrawWidget::mouseMoveEvent(QMouseEvent* pe)   // методы обработки события мыши при перемещении мыши
{
    if(catched_point == Point)
    {
        int dx = pe->x() - pr_mouseX;
        int dy = pe->y() - pr_mouseY;

        pr_mouseX = pe->position().x();
        pr_mouseY = pe->position().y();

        vector2d dd(dx, dy);
        for (int i = 0; i < int(point.size()); i++)
        {
            if (point[i] == k)
            {
                this->point[i] +=dd;
                k = point[i];
            }

        }

     this->repaint();
    }
}

void UDrawWidget::mouseReleaseEvent(QMouseEvent* pe) // методы обработки событий мыши при отжатии клавиши мыши
{
    catched_point = NoPoint;
}
