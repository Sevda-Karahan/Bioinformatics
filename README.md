1. Kullanılan Veri Seti

GSE19188

(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19188)

Bu proje, Homo sapiens organizmasında erken evre küçük hücreli olmayan akciğer
kanseri (NSCLC) ile ilgili verilerini kapsamaktadır. Veriler, gen ekspresyon analizi ile elde
edilmiştir.

Organizma: Homo sapiens

Deney Tipi: Expression profiling by array

Örnekler:
- 91 tümör örneği
- 65 bitişik normal akciğer doku örneği
- 
Bu veri seti, NSCLC'nin erken evresinde gen ekspresyon profillerini analiz ederek, tümör
ve histolojik imzaları tanımlamak amacıyla kullanılmıştır. Çalışma, farklı gen setlerinin
hastalığın teşhisi üzerindeki etkilerini belirlemeyi amaçlamaktadır.

2. Veri Setinin Hazırlanması ve R ile Ön İşleme

"GSE19188" veri setinin indirilmesi, işlenmesi, normalizasyonu gibi işlemleri “R” dilinde
gerçekleştirip hazırladığımız veri setini daha sonra kullanılmak üzere CSV formatında
kaydettik. 

3. Python ile Model Eğitimi

Kodlama ortamı olarak “Google Colab” kullanılmıştır

Ön işleme adımından elde edilen "GSE19188.csv" veri seti dosyası üzerinde SVM, KNN, Random Forest modelleri eğitilmiştir.
